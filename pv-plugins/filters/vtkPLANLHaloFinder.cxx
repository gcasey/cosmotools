/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkPLANLHaloFinder.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/*=========================================================================

  Program:   VTK/ParaView Los Alamos National Laboratory Modules (PVLANL)
  Module:    vtkPLANLHaloFinder.cxx

Copyright (c) 2007, 2009, Los Alamos National Security, LLC

All rights reserved.

Copyright 2007, 2009. Los Alamos National Security, LLC.
This software was produced under U.S. Government contract DE-AC52-06NA25396
for Los Alamos National Laboratory (LANL), which is operated by
Los Alamos National Security, LLC for the U.S. Department of Energy.
The U.S. Government has rights to use, reproduce, and distribute this software.
NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,
EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS SOFTWARE.
If software is modified to produce derivative works, such modified software
should be clearly marked, so as not to confuse it with the version available
from LANL.

Additionally, redistribution and use in source and binary forms, with or
without modification, are permitted provided that the following conditions
are met:
-   Redistributions of source code must retain the above copyright notice,
    this list of conditions and the following disclaimer.
-   Redistributions in binary form must reproduce the above copyright notice,
    this list of conditions and the following disclaimer in the documentation
    and/or other materials provided with the distribution.
-   Neither the name of Los Alamos National Security, LLC, Los Alamos National
    Laboratory, LANL, the U.S. Government, nor the names of its contributors
    may be used to endorse or promote products derived from this software
    without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY LOS ALAMOS NATIONAL SECURITY, LLC AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LOS ALAMOS NATIONAL SECURITY, LLC OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

=========================================================================*/

#ifndef USE_VTK_COSMO
#define USE_VTK_COSMO
#endif

#include "vtkPLANLHaloFinder.h"

// CosmologyTools includes
#include "CosmologyToolsMacros.h"

// VTK includes
#include "vtkDemandDrivenPipeline.h"
#include "vtkDummyController.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnstructuredGrid.h"

// halofinder includes
#include "ChainingMesh.h"
#include "CosmoDefinition.h"
#include "CosmoHaloFinderP.h"
#include "FOFHaloProperties.h"
#include "HaloCenterFinder.h"
#include "Partition.h"
#include "SODHalo.h"

// C++ includes
#include <cassert>
#include <vector>


vtkStandardNewMacro(vtkPLANLHaloFinder);

//------------------------------------------------------------------------------
vtkPLANLHaloFinder::vtkPLANLHaloFinder()
{
  this->SetNumberOfOutputPorts(2);

  this->Controller = 0;
  this->SetController(vtkMultiProcessController::GetGlobalController());

  this->NP = 256;
  this->RL = 100;
  this->Overlap = 5;
  this->BB = .2;
  this->PMin = 100;
  this->CopyHaloDataToParticles = 0;

  this->ComputeMostBoundParticle = 0;
  this->ComputeMostConnectedParticle = 0;

  this->ComputeSOD = 0;
  this->SODCenterType = 0;

  this->RhoC = RHO_C;
  this->SODMass = SOD_MASS;
  this->MinRadiusFactor = MIN_RADIUS_FACTOR;
  this->MaxRadiusFactor = MAX_RADIUS_FACTOR;
  this->SODBins = NUM_SOD_BINS;
  this->MinFOFSize = MIN_SOD_SIZE;
  this->MinFOFMass = MIN_SOD_MASS;

  this->HaloFinder = new CosmoHaloFinderP();
}

//------------------------------------------------------------------------------
vtkPLANLHaloFinder::~vtkPLANLHaloFinder()
{
  this->SetController( NULL );

  if( this->HaloFinder != NULL )
    {
    delete this->HaloFinder;
    }

  this->xx.clear();
  this->yy.clear();
  this->zz.clear();
  this->vx.clear();
  this->vy.clear();
  this->vz.clear();
  this->mass.clear();
  this->tag.clear();
  this->status.clear();
  this->potential.clear();
  this->mask.clear();
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if (this->Controller)
    {
    os << indent << "Controller: " << this->Controller << endl;
    }
  else
    {
    os << indent << "Controller: (null)\n";
    }
  os << indent << "NP: " << this->NP << endl;
  os << indent << "rL: " << this->RL << endl;
  os << indent << "Overlap: " << this->Overlap << endl;
  os << indent << "bb: " << this->BB << endl;
  os << indent << "pmin: " << this->PMin << endl;
  os << indent << "CopyHaloDataToParticles: "
     << this->CopyHaloDataToParticles << endl;
  os << indent << "ComputeMostBoundParticle: "
     << this->ComputeMostBoundParticle << endl;
  os << indent << "ComputeMostConnectedParticle: "
     << this->ComputeMostConnectedParticle << endl;
  os << indent << "ComputeSOD: "
     << this->ComputeSOD << endl;
  os << indent << "SODCenterType: "
     << this->SODCenterType << endl;
  os << indent << "RhoC: " << this->RhoC << endl;
  os << indent << "SODMass: " << this->SODMass << endl;
  os << indent << "MinRadiusFactor: " << this->MinRadiusFactor << endl;
  os << indent << "MaxRadiusFactor: " << this->MaxRadiusFactor << endl;
  os << indent << "SODBins: " << this->SODBins << endl;
  os << indent << "MinFOFSize: " << this->MinFOFSize << endl;
  os << indent << "MinFOFMass: " << this->MinFOFMass << endl;
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::SetController(vtkMultiProcessController *c)
{
  if(this->Controller == c)
    {
    return;
    }

  this->Modified();

  if(this->Controller != NULL)
    {
    this->Controller->UnRegister(this);
    this->Controller = 0;
    }

  if(c == NULL)
    {
    return;
    }

  this->Controller = c;
  c->Register(this);
}

//------------------------------------------------------------------------------
vtkMultiProcessController* vtkPLANLHaloFinder::GetController()
{
  return (vtkMultiProcessController*)this->Controller;
}

//------------------------------------------------------------------------------
int vtkPLANLHaloFinder::RequestInformation
(vtkInformation* vtkNotUsed(request),
 vtkInformationVector** inputVector,
 vtkInformationVector* outputVector)
{
  assert("pre: controller should not be NULL!" && (this->Controller != NULL));

  // set the other outputs to have the same number of pieces
  if((*inputVector)->GetInformationObject(0)->Has(
        vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()))
    {
    if(outputVector->GetInformationObject(1)->Has(
          vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()))
      {
      if(outputVector->GetInformationObject(0)->Get
         (vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()) !=
         outputVector->GetInformationObject(1)->Get
         (vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()))
        {
        outputVector->GetInformationObject(1)->Set
          (vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
           outputVector->GetInformationObject(0)->Get
           (vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()));
        }
      }
    else
      {
      outputVector->GetInformationObject(1)->Set
        (vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES(),
         outputVector->GetInformationObject(0)->Get
         (vtkStreamingDemandDrivenPipeline::MAXIMUM_NUMBER_OF_PIECES()));
      }
    }
  return 1;
}

//------------------------------------------------------------------------------
int vtkPLANLHaloFinder::RequestData(
          vtkInformation* request,
          vtkInformationVector** inputVector,
          vtkInformationVector* outputVector)
{
  assert("pre: controller should not be NULL!" && (this->Controller != NULL));

  // STEP 0: Get input object
  vtkInformation *input = inputVector[0]->GetInformationObject(0);
  assert( "pre: input information object is NULL " && (input != NULL) );
  vtkUnstructuredGrid *inputParticles =
      vtkUnstructuredGrid::SafeDownCast(
          input->Get(vtkDataObject::DATA_OBJECT() ) );
  assert("pre: input particles is NULL!" && (inputParticles != NULL) );

  // STEP 1: Get output objects. The output consists of two objects: (1) The
  // particles with halo information attached to it and (2) the halo centers
  // and generic FOF information.
  vtkInformation *output0 = outputVector->GetInformationObject(0);
  assert("pre: output information object is NULL" && (output0 != NULL) );
  vtkUnstructuredGrid *outputParticles =
      vtkUnstructuredGrid::SafeDownCast(
          output0->Get(vtkDataObject::DATA_OBJECT() ) );
  assert("pre: output particles is NULL" && (outputParticles != NULL) );

  vtkInformation *output1 = outputVector->GetInformationObject(1);
  assert("pre: output information object is NULL" && (output1 != NULL) );
  vtkUnstructuredGrid *haloCenters =
      vtkUnstructuredGrid::SafeDownCast(
          output1->Get(vtkDataObject::DATA_OBJECT() ) );
  assert("pre: halocenters is NULL" && (haloCenters != NULL) );

  // STEP 2: Get 1st output as a shallow-copy of the input & ensure integrity
  outputParticles->ShallowCopy( inputParticles );
  if( !this->CheckOutputIntegrity(outputParticles) )
    {
    vtkErrorMacro("Missing arrays from output particles mesh!");
    return 0;
    }

  // STEP 3: Initialize the partitioner used by the halo-finder which uses
  // MPI cartesian topology. Currently, the LANL halofinder assumes
  // MPI_COMM_WORLD. This should be changed in the short future.
  Partition::initialize();

  // STEP 4: Compute the FOF halos
  this->ComputeFOFHalos(outputParticles, haloCenters);

  // STEP 5: Compute SOD halos
  if( this->ComputeSOD )
    {
    // TODO: implement this
    }

  // STEP 6: Synchronize processes
  this->Controller->Barrier();
  return 1;
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::VectorizeData(
    vtkUnstructuredGrid *particles)
{
  assert("pre: input particles mesh is NULL" && (particles != NULL) );

  // TODO: VTK data arrays must match what POSVEL_T & ID_T are.

  vtkPoints     *points   = particles->GetPoints();
  vtkFloatArray *velocity = vtkFloatArray::SafeDownCast(
      particles->GetPointData()->GetArray("velocity") );
  vtkFloatArray *pmass = vtkFloatArray::SafeDownCast(
      particles->GetPointData()->GetArray("mass") );
  vtkIntArray *uid = vtkIntArray::SafeDownCast(
      particles->GetPointData()->GetArray("tag")  );
  vtkIntArray *owner = vtkIntArray::SafeDownCast(
      particles->GetPointData()->GetArray("ghost") );

  vtkIdType numParticles = points->GetNumberOfPoints();
  this->xx.resize( numParticles );
  this->yy.resize( numParticles );
  this->zz.resize( numParticles );
  this->vx.resize( numParticles );
  this->vy.resize( numParticles );
  this->vz.resize( numParticles );
  this->mass.resize( numParticles );
  this->tag.resize( numParticles );
  this->status.resize( numParticles );

  for( vtkIdType idx=0; idx < numParticles; ++idx )
    {
    // Extract position vector
    this->xx[ idx ] = points->GetPoint( idx )[ 0 ];
    this->yy[ idx ] = points->GetPoint( idx )[ 1 ];
    this->zz[ idx ] = points->GetPoint( idx )[ 2 ];

    // Extract velocity vector
    this->vx[ idx ] = velocity->GetComponent(idx, 0);
    this->vy[ idx ] = velocity->GetComponent(idx, 1);
    this->vz[ idx ] = velocity->GetComponent(idx, 2);

    // Extract global particle ID information
    this->tag[ idx ] = uid->GetValue( idx );

    // Extract status
    this->status[ idx ] = owner->GetValue( idx );
    } // END for all particles

  // Remove the ghost array since the status vector now has that information
  particles->GetPointData()->RemoveArray("ghost");
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::ComputeFOFHalos(
    vtkUnstructuredGrid *particles, vtkUnstructuredGrid *haloCenters)
{
  assert("pre: input particles mesh is NULL" && (particles != NULL) );
  assert("pre: halo-centers data-structure is NULL" && (haloCenters != NULL) );
  assert("pre: halo-finder is not allocated!" && (this->HaloFinder != NULL) );

  // STEP 0: Vectorize the data
  this->VectorizeData(particles);

  // STEP 1: Allocate potential & mask vectors required by the halo-finder
  this->potential.resize( particles->GetNumberOfPoints() );
  this->mask.resize( particles->GetNumberOfPoints() );

  // STEP 2: Initialize halo-finder parameters
  this->HaloFinder->setParameters(
      "",this->RL,this->Overlap,this->NP,this->PMin,this->BB);
  this->HaloFinder->setParticles(
      &this->xx,&this->yy,&this->zz,
      &this->vx,&this->vy,&this->vz,
      &this->potential,
      &this->tag,
      &this->mask,
      &this->status
      );

  // STEP 3: Execute the halo-finder
  this->HaloFinder->executeHaloFinder();
  this->HaloFinder->collectHalos();
  this->HaloFinder->mergeHalos();

  // STEP 4: Calculate basic FOF halo properties
  std::vector<POSVEL_T> fofMass; /* mass of every halo */
  std::vector<POSVEL_T> fofXPos; /* x-component of the FOF position */
  std::vector<POSVEL_T> fofYPos; /* y-component of the FOF position */
  std::vector<POSVEL_T> fofZPos; /* z-component of the FOF position */
  std::vector<POSVEL_T> fofXVel; /* x-component of the FOF velocity */
  std::vector<POSVEL_T> fofYVel; /* y-component of the FOF velocity */
  std::vector<POSVEL_T> fofZVel; /* z-component of the FOF velocity */
  std::vector<POSVEL_T> fofXCofMass;/* x-component of the FOF center of mass */
  std::vector<POSVEL_T> fofYCofMass;/* y-component of the FOF center of mass */
  std::vector<POSVEL_T> fofZCofMass;/* z-component of the FOF center of mass */
  std::vector<POSVEL_T> fofVelDisp; /* velocity dispersion of every halo */
  this->ComputeFOFHaloProperties(
      fofMass,
      fofXPos,fofYPos,fofZPos,
      fofXVel,fofYVel,fofZVel,
      fofXCofMass,fofYCofMass,fofZCofMass,
      fofVelDisp
      );



}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::ComputeFOFHaloProperties(
      std::vector<POSVEL_T> &fofMass,
      std::vector<POSVEL_T> &fofXPos,
      std::vector<POSVEL_T> &fofYPos,
      std::vector<POSVEL_T> &fofZPos,
      std::vector<POSVEL_T> &fofXVel,
      std::vector<POSVEL_T> &fofYVel,
      std::vector<POSVEL_T> &fofZVel,
      std::vector<POSVEL_T> &fofXCofMass,
      std::vector<POSVEL_T> &fofYCofMass,
      std::vector<POSVEL_T> &fofZCofMass,
      std::vector<POSVEL_T> &fofVelDisp)
{
  assert("pre: haloFinder is NULL!" && (this->HaloFinder != NULL) );

  int numberOfHalos = this->HaloFinder->getNumberOfHalos();
  int* fofHalos     = this->HaloFinder->getHalos();
  int* fofHaloCount = this->HaloFinder->getHaloCount();
  int* fofHaloList  = this->HaloFinder->getHaloList();

  FOFHaloProperties* fof = new FOFHaloProperties();
  fof->setHalos(numberOfHalos,fofHalos,fofHaloCount,fofHaloList);
  fof->setParameters("",this->RL,this->Overlap,this->BB);
  fof->setParticles(
      &this->xx,&this->yy,&this->zz,
      &this->vx,&this->vy,&this->vz,
      &this->mass,
      &this->potential,
      &this->tag,
      &this->mask,
      &this->status);

  // Compute average halo position
  fof->FOFPosition(&fofXPos,&fofYPos,&fofZPos);

  // Compute mass of every halo
  fof->FOFHaloMass(&fofMass);

  // Compute center of mass of every FOF halo
  fof->FOFCenterOfMass(&fofXCofMass,&fofYCofMass,&fofZCofMass);

  // Compute average velocity of every FOF halo
  fof->FOFVelocity(&fofXVel,&fofYVel,&fofZVel);

  // Compute the velocity dispersion
  fof->FOFVelocityDispersion(&fofXVel,&fofYVel,&fofZVel,&fofVelDisp);
  delete fof;
}

//------------------------------------------------------------------------------
bool vtkPLANLHaloFinder::CheckOutputIntegrity(
      vtkUnstructuredGrid* outputParticles)
{
  assert("pre: particle mesh is NULL" && (outputParticles != NULL) );
  if(!outputParticles->GetPointData()->HasArray("velocity") ||
      !outputParticles->GetPointData()->HasArray("mass") ||
      !outputParticles->GetPointData()->HasArray("tag") ||
      !outputParticles->GetPointData()->HasArray("ghost"))
      {
      vtkErrorMacro(
          << "The input data does not have one or more of "
          << "the following point arrays: velocity, mass, tag, or ghost."
         );
      return false;
      }
  return true;
}


