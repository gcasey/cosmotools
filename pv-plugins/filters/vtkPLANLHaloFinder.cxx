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
#include "vtkCellType.h"
#include "vtkDemandDrivenPipeline.h"
#include "vtkDoubleArray.h"
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

  this->NP      = 256;
  this->RL      = 100;
  this->Overlap = 5;
  this->BB      = .2;
  this->PMin    = 100;

  this->ComputeSOD          = 0;
  this->CenterFindingMethod = AVERAGE;

  this->HaloFinder      = NULL;
  this->RhoC            = cosmologytools::RHO_C;
  this->SODMass         = cosmologytools::SOD_MASS;
  this->MinRadiusFactor = cosmologytools::MIN_RADIUS_FACTOR;
  this->MaxRadiusFactor = cosmologytools::MAX_RADIUS_FACTOR;
  this->SODBins         = cosmologytools::NUM_SOD_BINS;
  this->MinFOFSize      = cosmologytools::MIN_SOD_SIZE;
  this->MinFOFMass      = cosmologytools::MIN_SOD_MASS;
}

//------------------------------------------------------------------------------
vtkPLANLHaloFinder::~vtkPLANLHaloFinder()
{
  this->Controller = NULL;

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

  this->fofMass.clear();
  this->fofXPos.clear();
  this->fofYPos.clear();
  this->fofZPos.clear();
  this->fofXVel.clear();
  this->fofYVel.clear();
  this->fofZVel.clear();
  this->fofXCofMass.clear();
  this->fofYCofMass.clear();
  this->fofZCofMass.clear();
  this->fofVelDisp.clear();

  this->ExtractedHalos.clear();
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::SetController(vtkMultiProcessController *c)
{
  assert("pre: cannot set a NULL controller!" && (c != NULL) );
  this->Controller = c;
}

//------------------------------------------------------------------------------
vtkMultiProcessController* vtkPLANLHaloFinder::GetController()
{
  return(this->Controller);
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

  // Reset previously calculated data
  this->ResetHaloFinderInternals();

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

  if( inputParticles->GetNumberOfPoints() == 0 )
    {
    // Empty input
    return 1;
    }

  // STEP 2: Get 1st output as a shallow-copy of the input & ensure integrity
  outputParticles->ShallowCopy( inputParticles );
  if( !this->CheckOutputIntegrity(outputParticles) )
    {
    vtkErrorMacro("Missing arrays from output particles mesh!");
    return 0;
    }


  // STEP 3: Initialize the partitioner used by the halo-finder which uses
  // MPI cartesian topology. Currently, the LANL halofinder assumes
  // MPI_COMM_WORLD!!!!!!. This should be changed in the short future.
  cosmologytools::Partition::initialize();

  // Delete previously computed results to re-compute with modified params
  // TODO: We need to get smarter about this in the future, e.g., determine
  // if we really have to run the halo-finder again or just filter differently
  // the results, just change the halo-centers etc.
  if( this->HaloFinder != NULL )
    {
    delete this->HaloFinder;
    }
  this->HaloFinder = new cosmologytools::CosmoHaloFinderP();

  // STEP 4: Compute the FOF halos
  this->ComputeFOFHalos(outputParticles, haloCenters);

  // STEP 5: Compute SOD halos
  if( this->ComputeSOD )
    {
    this->ComputeSODHalos(outputParticles, haloCenters);
    }

  // STEP 6: Synchronize processes
  this->Controller->Barrier();
  return 1;
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::ComputeSODHalos(
      vtkUnstructuredGrid *particles,
      vtkUnstructuredGrid *fofHaloCenters)
{
  assert("pre: input particles should not be NULL" &&
         (particles != NULL) );
  assert("pre: FOF halo-centers should not be NULL" &&
         (fofHaloCenters != NULL));

  // STEP 0: Initialize SOD arrays & acquire handles
  // "SODAveragePosition"
  // "SODCenterOfMass"
  // "SODMass"
  // "SODAverageVelocity"
  // "SODVelocityDispersion"
  // "SODRadius"
  this->InitializeSODHaloArrays( fofHaloCenters );
  vtkPointData *PD = fofHaloCenters->GetPointData();
  assert("pre: missing SODAveragePosition array" &&
          PD->HasArray("SODAveragePosition") );
  assert("pre: missing SODCenterOfMass" &&
          PD->HasArray("SODCenterOfMass") );
  assert("pre: missing SODMass" &&
          PD->HasArray("SODMass") );
  assert("pre: missing SODAverageVelocity" &&
          PD->HasArray("SODAverageVelocity") );
  assert("pre: missing SODVelocityDispersion" &&
          PD->HasArray("SODVelocityDispersion") );
  assert("pre: missing SODRadius" &&
          PD->HasArray("SODRadius") );
  vtkDoubleArray *sodPos =
      vtkDoubleArray::SafeDownCast(PD->GetArray("SODAveragePosition"));
  vtkDoubleArray *sodCofMass =
      vtkDoubleArray::SafeDownCast(PD->GetArray("SODCenterOfMass"));
  vtkDoubleArray *sodMass =
      vtkDoubleArray::SafeDownCast(PD->GetArray("SODMass"));
  vtkDoubleArray *sodVelocity =
      vtkDoubleArray::SafeDownCast(PD->GetArray("SODAverageVelocity"));
  vtkDoubleArray *sodDispersion =
      vtkDoubleArray::SafeDownCast(PD->GetArray("SODVelocityDispersion"));
  vtkDoubleArray *sodRadius =
      vtkDoubleArray::SafeDownCast(PD->GetArray("SODRadius"));

  // STEP 1: Construct the ChainingMesh
  cosmologytools::ChainingMesh *chainMesh =
      new cosmologytools::ChainingMesh(
          this->RL,this->Overlap,cosmologytools::CHAIN_SIZE,
          &this->xx,&this->yy,&this->zz);

  // STEP 2: Loop through all halos and compute SOD halos
  for(unsigned int i=0; i < this->ExtractedHalos.size(); ++i)
    {
    int internalHaloIdx = this->ExtractedHalos[ i ];
    int haloSize        = this->HaloFinder->getHaloCount()[internalHaloIdx];

    double haloMass     = this->fofMass[internalHaloIdx];

    if( (haloMass < this->MinFOFMass) || (haloSize < this->MinFOFSize) )
      {
      continue;
      }

    cosmologytools::SODHalo *sod = new cosmologytools::SODHalo();
    sod->setParameters( chainMesh, this->SODBins, this->RL, this->NP,
                        this->RhoC, this->SODMass, this->RhoC,
                        this->MinRadiusFactor, this->MaxRadiusFactor );
    sod->setParticles(
        &this->xx,&this->yy,&this->zz,
        &this->vx,&this->vy,&this->vz,
        &this->mass,
        &this->tag);

    double center[3];
    fofHaloCenters->GetPoint(i,center);
    sod->createSODHalo(
        this->HaloFinder->getHaloCount()[internalHaloIdx],
        center[0],center[1],center[2],
        this->fofXVel[internalHaloIdx],
        this->fofYVel[internalHaloIdx],
        this->fofZVel[internalHaloIdx],
        this->fofMass[internalHaloIdx]);

    if(sod->SODHaloSize() > 0)
      {
      POSVEL_T pos[3];
      POSVEL_T cofmass[3];
      POSVEL_T mass;
      POSVEL_T vel[3];
      POSVEL_T disp;
      POSVEL_T radius = sod->SODRadius();

      sod->SODAverageLocation(pos);
      sod->SODCenterOfMass(cofmass);
      sod->SODMass(&mass);
      sod->SODAverageVelocity(vel);
      sod->SODVelocityDispersion(&disp);

      sodPos->SetTuple3(i,pos[0],pos[1],pos[2]);
      sodCofMass->SetTuple3(i,cofmass[0],cofmass[1],cofmass[2]);
      sodMass->SetValue(i,mass);
      sodVelocity->SetTuple3(i,vel[0],vel[1],vel[2]);
      sodDispersion->SetValue(i,disp);
      sodRadius->SetValue(i,radius);
      }

    delete sod;
    } // END for all halos within the PMIN threshold

  // STEP 3: De-allocate Chain mesh
  delete chainMesh;
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::InitializeSODHaloArrays(
      vtkUnstructuredGrid *haloCenters)
{
  assert("pre: centers is NULL" && (haloCenters != NULL) );

  //TODO: Again, these arrays should match the type of POSVEL_T

  // STEP 0: Allocate arrays
  vtkDoubleArray *averagePosition = vtkDoubleArray::New();
  averagePosition->SetName("SODAveragePosition");
  averagePosition->SetNumberOfComponents( 3 );
  averagePosition->SetNumberOfTuples( haloCenters->GetNumberOfPoints() );

  vtkDoubleArray *centerOfMass = vtkDoubleArray::New();
  centerOfMass->SetName("SODCenterOfMass");
  centerOfMass->SetNumberOfComponents( 3 );
  centerOfMass->SetNumberOfTuples( haloCenters->GetNumberOfPoints() );

  vtkDoubleArray *mass = vtkDoubleArray::New();
  mass->SetName("SODMass");
  mass->SetNumberOfComponents( 1 );
  mass->SetNumberOfTuples( haloCenters->GetNumberOfPoints() );

  vtkDoubleArray *averageVelocity = vtkDoubleArray::New();
  averageVelocity->SetName("SODAverageVelocity");
  averageVelocity->SetNumberOfComponents( 3 );
  averageVelocity->SetNumberOfTuples( haloCenters->GetNumberOfPoints() );

  vtkDoubleArray *velDispersion = vtkDoubleArray::New();
  velDispersion->SetName("SODVelocityDispersion");
  velDispersion->SetNumberOfComponents( 1 );
  velDispersion->SetNumberOfTuples( haloCenters->GetNumberOfPoints() );

  vtkDoubleArray *radius = vtkDoubleArray::New();
  radius->SetName("SODRadius");
  radius->SetNumberOfComponents( 1 );
  radius->SetNumberOfTuples( haloCenters->GetNumberOfPoints() );

  // STEP 1: Initialize
  for( vtkIdType halo=0; halo < haloCenters->GetNumberOfPoints(); ++halo )
    {
    averagePosition->SetTuple3(halo,0.0,0.0,0.0);
    centerOfMass->SetTuple3(halo,0.0,0.0,0.0);
    mass->SetValue(halo,0.0);
    averageVelocity->SetTuple3(halo,0.0,0.0,0.0);
    velDispersion->SetValue(halo,0.0);
    radius->SetValue(halo,0.0);
    } // END for all extracted FOF halos

  // STEP 2: Add arrays to halo-centers
  haloCenters->GetPointData()->AddArray(averagePosition);
  averagePosition->Delete();
  haloCenters->GetPointData()->AddArray(centerOfMass);
  centerOfMass->Delete();
  haloCenters->GetPointData()->AddArray(mass);
  mass->Delete();
  haloCenters->GetPointData()->AddArray(averageVelocity);
  averageVelocity->Delete();
  haloCenters->GetPointData()->AddArray(velDispersion);
  velDispersion->Delete();
  haloCenters->GetPointData()->AddArray(radius);
  radius->Delete();
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

  vtkIntArray *haloTag = vtkIntArray::New();
  haloTag->SetName( "HaloID" );
  haloTag->SetNumberOfComponents( 1 );
  haloTag->SetNumberOfTuples( points->GetNumberOfPoints( ) );
  int *haloTagPtr = static_cast<int*>( haloTag->GetVoidPointer(0) );

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

    // Extract the mass
    this->mass[ idx ] = pmass->GetValue( idx );

    // Extract global particle ID information & also setup global-to-local map
    this->tag[ idx ] = uid->GetValue( idx );

    // Extract status
    this->status[ idx ] = owner->GetValue( idx );

    // Initialize all halo IDs to -1, i.e., all particles are not in halos
    haloTagPtr[ idx ] = -1;
    } // END for all particles

  // Remove the ghost array since the status vector now has that information
  particles->GetPointData()->RemoveArray("ghost");

  // Add halo IDs to the particles
  particles->GetPointData()->AddArray( haloTag );
  haloTag->Delete();
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
  this->ComputeFOFHaloProperties();

  // STEP 5: Filter out halos within the PMin threshold
  int numberOfFOFHalos = this->HaloFinder->getNumberOfHalos();
  int *fofHaloCount    = this->HaloFinder->getHaloCount();

  for( int halo=0; halo < numberOfFOFHalos; ++halo )
    {
    // Disregard halos that do not fall within the PMin threshold
    if( fofHaloCount[ halo ] >= this->PMin )
      {
      this->ExtractedHalos.push_back( halo );
      } // END if haloSize is within threshold
    } // END for all halos

  // STEP 6: Loop through the extracted halos and do the following:
  //          1. Compute the halo-centers
  //          2. Attach halo-properties to each halo-center,e.g.,mass,vel,etc.
  //          3. Mark all particles within each halo
  this->InitializeHaloCenters( haloCenters, this->ExtractedHalos.size() );
  assert("pre: halo mass array not present!" &&
          haloCenters->GetPointData()->HasArray("HaloMass") );
  assert("pre: AverageVelocity array not present!" &&
          haloCenters->GetPointData()->HasArray("AverageVelocity") );
  assert("pre: VelocityDispersion array not present!" &&
          haloCenters->GetPointData()->HasArray("VelocityDispersion") );
  assert("pre: HaloID array not present!" &&
          haloCenters->GetPointData()->HasArray("HaloID") );

  vtkPoints *pnts  = haloCenters->GetPoints();
  vtkPointData *PD = haloCenters->GetPointData();
  double* haloMass = static_cast<double*>(
      PD->GetArray("HaloMass")->GetVoidPointer(0));
  double* haloAverageVel = static_cast<double*>(
      PD->GetArray("AverageVelocity")->GetVoidPointer(0));
  double* haloVelDisp = static_cast<double*>(
      PD->GetArray("VelocityDispersion")->GetVoidPointer(0));
  int* haloId = static_cast<int*>(
      PD->GetArray("HaloID")->GetVoidPointer(0));

  double center[3];
  for(unsigned int halo=0; halo < this->ExtractedHalos.size(); ++halo )
    {
    int haloIdx = this->ExtractedHalos[ halo ];
    assert( "pre: haloIdx is out-of-bounds!" &&
            (haloIdx >= 0) && (haloIdx < this->fofMass.size() ) );

    this->MarkHaloParticlesAndGetCenter(halo,haloIdx,center,particles);
    pnts->SetPoint( halo, center );

    haloMass[ halo ]          = this->fofMass[ haloIdx ];
    haloVelDisp[ halo ]       = this->fofVelDisp[ haloIdx ];
    haloAverageVel[ halo*3  ] = this->fofXVel[ haloIdx ];
    haloAverageVel[ halo*3+1] = this->fofYVel[ haloIdx ];
    haloAverageVel[ halo*3+2] = this->fofZVel[ haloIdx ];
    haloId[ halo ]            = halo;
    } // END for all extracted halos
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::MarkHaloParticlesAndGetCenter(
      const unsigned int halo,
      const int internalHaloIdx,
      double center[3],
      vtkUnstructuredGrid *particles)
{
  assert("pre: output particles data-structured is NULL" &&
         (particles != NULL) );
  assert("pre: particles must have a HaloID array" &&
         (particles->GetPointData()->HasArray("HaloID")));

  // STEP 0: Get pointer to the haloTags array
  int *haloTags = static_cast<int*>(
      particles->GetPointData()->GetArray("HaloID")->GetVoidPointer(0));

  // STEP 1: Construct FOF halo properties (TODO: we can modularize this part)
  int numberOfHalos = this->HaloFinder->getNumberOfHalos();
  int* fofHalos     = this->HaloFinder->getHalos();
  int* fofHaloCount = this->HaloFinder->getHaloCount();
  int* fofHaloList  = this->HaloFinder->getHaloList();

  cosmologytools::FOFHaloProperties* fof =
      new cosmologytools::FOFHaloProperties();
  fof->setHalos(numberOfHalos,fofHalos,fofHaloCount,fofHaloList);
  fof->setParameters("",this->RL,this->Overlap,this->BB);
  fof->setParticles(
      &this->xx,&this->yy,&this->zz,
      &this->vx,&this->vy,&this->vz,
      &this->mass,
      &this->potential,
      &this->tag,
      &this->mask,
      &this->status
      );

  // STEP 2: Get the particle halo information for the given halo with the given
  // internal halo idx. The halo finder find halos in [0, N]. However, we
  // filter out the halos that are not within the PMin threshold. This yields
  // halos [0, M] where M < N. The internalHaloIdx is w.r.t. the [0,N] range
  // while the halo is w.r.t. the [0,M] range.
  long size = fofHaloCount[ internalHaloIdx ];
  POSVEL_T* xLocHalo = new POSVEL_T[size];
  POSVEL_T* yLocHalo = new POSVEL_T[size];
  POSVEL_T* zLocHalo = new POSVEL_T[size];
  POSVEL_T* xVelHalo = new POSVEL_T[size];
  POSVEL_T* yVelHalo = new POSVEL_T[size];
  POSVEL_T* zVelHalo = new POSVEL_T[size];
  POSVEL_T* massHalo = new POSVEL_T[size];
  ID_T* id           = new ID_T[size];
  int* actualIdx = new int[size];
  fof->extractInformation(
      internalHaloIdx,actualIdx,
      xLocHalo, yLocHalo, zLocHalo,
      xVelHalo, yVelHalo, zVelHalo,
      massHalo, id
      );

  // STEP 3: Mark particles within the halo with the the given halo tag/ID.
  for( int haloParticleIdx=0; haloParticleIdx < size; ++haloParticleIdx )
    {
    haloTags[ actualIdx[ haloParticleIdx ] ] = halo;
    } // END for all haloParticles

  // STEP 4: Find the center
  switch( this->CenterFindingMethod )
    {
    case CENTER_OF_MASS:
      assert("pre:center of mass has not been constructed correctly!" &&
             (this->fofXCofMass.size()==numberOfHalos));
      assert("pre:center of mass has not been constructed correctly!" &&
             (this->fofYCofMass.size()==numberOfHalos));
      assert("pre:center of mass has not been constructed correctly!" &&
             (this->fofZCofMass.size()==numberOfHalos));
      center[0] = this->fofXCofMass[ internalHaloIdx ];
      center[1] = this->fofYCofMass[ internalHaloIdx ];
      center[2] = this->fofZCofMass[ internalHaloIdx ];
      break;
    case AVERAGE:
      assert("pre:average halo center has not been constructed correctly!" &&
             (this->fofXPos.size()==numberOfHalos));
      assert("pre:average halo center has not been constructed correctly!" &&
             (this->fofYPos.size()==numberOfHalos));
      assert("pre:average halo center has not been constructed correctly!" &&
             (this->fofZPos.size()==numberOfHalos));
      center[0] = this->fofXPos[ internalHaloIdx ];
      center[1] = this->fofYPos[ internalHaloIdx ];
      center[2] = this->fofZPos[ internalHaloIdx ];
      break;
    case MBP:
    case MCP:
      {
      int centerIndex = -1;
      POTENTIAL_T minPotential;

      cosmologytools::HaloCenterFinder centerFinder;
      centerFinder.setParticles(size,xLocHalo,yLocHalo,zLocHalo,massHalo,id);
      centerFinder.setParameters(this->BB,this->Overlap);

      if( this->CenterFindingMethod == MBP )
        {
        centerIndex = (size < cosmologytools::MBP_THRESHOLD)?
            centerFinder.mostBoundParticleN2( &minPotential ):
              centerFinder.mostBoundParticleAStar( &minPotential );
        } // END if MBP
      else
        {
        centerIndex = (size < cosmologytools::MCP_THRESHOLD)?
            centerFinder.mostConnectedParticleN2():
              centerFinder.mostConnectedParticleChainMesh();
        }

      int pidx = actualIdx[ centerIndex ];
      assert("pre: particle index is out-of-bounds" &&
             (pidx >= 0) && (pidx < particles->GetNumberOfPoints() ) );
      particles->GetPoint( pidx, center );
      }
      break;
    default:
      vtkErrorMacro("Undefined center-finding method!");
    }

  delete [] xLocHalo;
  delete [] yLocHalo;
  delete [] zLocHalo;
  delete [] xVelHalo;
  delete [] yVelHalo;
  delete [] zVelHalo;
  delete [] massHalo;
  delete [] id;
  delete [] actualIdx;
  delete fof;
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::InitializeHaloCenters(
    vtkUnstructuredGrid *haloCenters, unsigned int N)
{
  // TODO: Note the vtkArrays here should match what POSVEL_T, ID_T are set
  haloCenters->Allocate( N, 0);

  vtkPoints* pnts = vtkPoints::New();
  pnts->SetDataTypeToDouble();
  pnts->SetNumberOfPoints( N );

  vtkDoubleArray *haloMass = vtkDoubleArray::New();
  haloMass->SetName( "HaloMass" );
  haloMass->SetNumberOfComponents( 1 );
  haloMass->SetNumberOfTuples( N );
  double *haloMassArray =
      static_cast<double*>(haloMass->GetVoidPointer(0));

  vtkDoubleArray *averageVelocity = vtkDoubleArray::New();
  averageVelocity->SetName("AverageVelocity");
  averageVelocity->SetNumberOfComponents( 3 );
  averageVelocity->SetNumberOfTuples( N );
  double *velArray =
      static_cast<double*>(averageVelocity->GetVoidPointer(0));

  vtkDoubleArray *velDispersion = vtkDoubleArray::New();
  velDispersion->SetName("VelocityDispersion");
  velDispersion->SetNumberOfComponents( 1 );
  velDispersion->SetNumberOfTuples( N );
  double *velDispArray =
      static_cast<double*>(velDispersion->GetVoidPointer(0));

  vtkIntArray *haloIdx = vtkIntArray::New();
  haloIdx->SetName( "HaloID" );
  haloIdx->SetNumberOfComponents( 1 );
  haloIdx->SetNumberOfTuples( N );
  int *haloIdxArray =
      static_cast< int* >(haloIdx->GetVoidPointer(0));

  for(unsigned int i=0; i < N; ++i )
    {
    vtkIdType vertexIdx = i;
    haloCenters->InsertNextCell(VTK_VERTEX,1,&vertexIdx);
    velArray[ i*3 ]    =
    velArray[ i*3+1 ]  =
    velArray[ i*3+2 ]  =
    velDispArray[ i ]  = 0.0;
    haloIdxArray[ i ]  = -1;
    haloMassArray[ i ] = 0.0;
    }

  haloCenters->SetPoints( pnts );
  pnts->Delete();
  haloCenters->GetPointData()->AddArray(averageVelocity);
  averageVelocity->Delete();
  haloCenters->GetPointData()->AddArray(velDispersion);
  velDispersion->Delete();
  haloCenters->GetPointData()->AddArray(haloIdx);
  haloIdx->Delete();
  haloCenters->GetPointData()->AddArray(haloMass);
  haloMass->Delete();
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::ComputeFOFHaloProperties()
{
  assert("pre: haloFinder is NULL!" && (this->HaloFinder != NULL) );

  int numberOfHalos = this->HaloFinder->getNumberOfHalos();
  int* fofHalos     = this->HaloFinder->getHalos();
  int* fofHaloCount = this->HaloFinder->getHaloCount();
  int* fofHaloList  = this->HaloFinder->getHaloList();

  cosmologytools::FOFHaloProperties* fof =
      new cosmologytools::FOFHaloProperties();
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

  // Compute average halo position if that's what will be used as the halo
  // center position
  if( this->CenterFindingMethod == AVERAGE )
    {
    fof->FOFPosition(&this->fofXPos,&this->fofYPos,&this->fofZPos);
    }
  // Compute center of mass of every FOF halo if that's what will be used as
  // the halo center
  else if( this->CenterFindingMethod == CENTER_OF_MASS )
    {
    fof->FOFCenterOfMass(
      &this->fofXCofMass,&this->fofYCofMass,&this->fofZCofMass);
    }

  fof->FOFAttributes(
      &this->fofMass,
      &this->fofXVel, &this->fofYVel,&this->fofZVel,
      &this->fofVelDisp);

  // Sanity checks!
  assert("post: FOF mass property not correctly computed!" &&
         (this->fofMass.size()==numberOfHalos) );
  assert("post: FOF x-velocity component not correctly computed!" &&
         (this->fofXVel.size()==numberOfHalos) );
  assert("post: FOF y-velocity component not correctly computed!" &&
         (this->fofYVel.size()==numberOfHalos) );
  assert("post: FOF z-velocity component not correctly computed!" &&
         (this->fofZVel.size()==numberOfHalos) );
  assert("post: FOF velocity dispersion not correctly computed!" &&
         (this->fofVelDisp.size()==numberOfHalos));

  delete fof;
}

//------------------------------------------------------------------------------
void vtkPLANLHaloFinder::ResetHaloFinderInternals()
{
  // input particle information
  this->xx.resize(0); this->yy.resize(0);  this->zz.resize(0);
  this->vx.resize(0); this->vy.resize(0);  this->vz.resize(0);
  this->mass.resize(0);
  this->potential.resize(0);
  this->tag.resize(0);
  this->status.resize(0);
  this->mask.resize(0);

  // computed FOF properties
  this->fofMass.resize(0);
  this->fofXPos.resize(0); this->fofYPos.resize(0); this->fofZPos.resize(0);
  this->fofXVel.resize(0); this->fofYVel.resize(0); this->fofZVel.resize(0);
  this->fofXCofMass.resize(0);
  this->fofYCofMass.resize(0);
  this->fofZCofMass.resize(0);
  this->fofVelDisp.resize(0);
  this->ExtractedHalos.resize(0);
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


