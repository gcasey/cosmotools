/*=========================================================================

 Program:   Visualization Toolkit
 Module:    vtkCosmoDensityProfile.cxx

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/
#include "vtkCosmoDensityProfile.h"
#include "vtkObjectFactory.h"
#include "vtkIndent.h"
#include "vtkMultiProcessController.h"
#include "vtkMultiBlockDataSet.h"
#include "vtkTable.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkUnstructuredGrid.h"
#include "vtkSphere.h"
#include "vtkSphereSource.h"
#include "vtkIdTypeArray.h"
#include "vtkCommunicator.h"

#include <set>
#include <list>
#include <cassert>
#include <sstream>
#include <fstream>

namespace {

class vtkVectorSumOperator : public vtkCommunicator::Operation
{
  // Description:
  // Performs a "B+=A" vector operation.
  virtual void Function(
      const void *A, void *B, vtkIdType length, int datatype )
    {
      assert("pre: A vector is NULL" && (A != NULL) );
      assert("pre: B vector is NULL");
      vtkIdType idx = 0;
      switch( datatype )
        {
        case VTK_INT:
          {
          const int *aPtr = reinterpret_cast<const int*>(A);
          int *bPtr       = reinterpret_cast<int*>(B);
          for( ; idx < length; ++idx )
            {
            bPtr[ idx ] += aPtr[ idx ];
            }
          }
          break;
        case VTK_DOUBLE:
          {
          const double *aPtr = reinterpret_cast<const double*>(A);
          double *bPtr       = reinterpret_cast<double*>(B);
          for( ; idx < length; ++idx )
            {
            bPtr[ idx ] += aPtr[ idx ];
            }
          }
          break;
        default:
         std::cerr << "ERROR: Unresolved type! ";
         std::cerr << "Only VTK_INT and VTK_DOUBLE are supported\n";
        }
    }

  // Description:
  // Sets Commutative to true for this operation
  virtual int Commutative() { return 1; }
};

}

vtkStandardNewMacro(vtkCosmoDensityProfile);

//------------------------------------------------------------------------------
vtkCosmoDensityProfile::vtkCosmoDensityProfile()
{
  this->Controller = vtkMultiProcessController::GetGlobalController();
  this->Center[0]  = this->Center[1] = this->Center[2] = 0.0;
  this->Radius     = 5.0;
  this->NumberOfSpheres = 0;
  this->SetNumberOfInputPorts( 1 );
  this->SetNumberOfOutputPorts( 2 );
}

//------------------------------------------------------------------------------
vtkCosmoDensityProfile::~vtkCosmoDensityProfile()
{
  this->ConcentricRadii.clear();
  this->NumParticlesInSphere.clear();
}

//------------------------------------------------------------------------------
void vtkCosmoDensityProfile::PrintSelf(ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

//------------------------------------------------------------------------------
int vtkCosmoDensityProfile::FillInputPortInformation(
                      int vtkNotUsed(port),vtkInformation *info)
{
  assert("pre: input information is NULL!" && (info != NULL) );
  info->Set(
      vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkUnstructuredGrid");
  return 1;
}

//------------------------------------------------------------------------------
int vtkCosmoDensityProfile::FillOutputPortInformation(
                      int port,vtkInformation *info)
{
  switch( port )
    {
    case 0:
      info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkMultiBlockDataSet");
      break;
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkTable");
      break;
    default:
      vtkErrorMacro("Undefined port!" << port );
    }
  return 1;
}

//------------------------------------------------------------------------------
void vtkCosmoDensityProfile::GenerateSpheres(vtkMultiBlockDataSet *mbds)
{
  assert("pre: output multi-block dataset is NULL" && (mbds != NULL) );

  mbds->SetNumberOfBlocks( this->NumberOfSpheres );

  this->NumParticlesInSphere.resize( this->NumberOfSpheres,0 );
  this->ConcentricRadii.resize( this->NumberOfSpheres );
  double h = static_cast<double>( this->Radius/this->NumberOfSpheres );


  for( int i=0; i < this->NumberOfSpheres; ++i )
    {
    this->ConcentricRadii[ i ] = (i+1)*h;

    vtkSphereSource *sphereGenerator = vtkSphereSource::New();
    sphereGenerator->SetCenter( this->Center );
    sphereGenerator->SetRadius( this->ConcentricRadii[ i ] );
    sphereGenerator->Update();
    mbds->SetBlock( i, sphereGenerator->GetOutput() );
    sphereGenerator->Delete();
    }

}

//------------------------------------------------------------------------------
void vtkCosmoDensityProfile::FindParticlesInSphere(
    const int sphereIdx,
    vtkUnstructuredGrid* particles,
    std::list< vtkIdType > &particleIds )
{
  double r = this->ConcentricRadii[ sphereIdx ];
  vtkSphere *mySphere = vtkSphere::New();
  mySphere->SetCenter( this->Center );
  mySphere->SetRadius( r );

  double pnt[3];
  if( sphereIdx == this->NumberOfSpheres-1 )
    {
    for( vtkIdType idx=0; idx < particles->GetNumberOfPoints(); ++idx )
      {
      particles->GetPoint( idx, pnt );
      if( mySphere->EvaluateFunction(pnt) <= 0.0 )
        {
        particleIds.push_front( idx );
        }
      } // END for all particles
    }
  else
    {
    std::set< vtkIdType > particlesToErase;
    std::list< vtkIdType >::iterator iter = particleIds.begin();
    for( ; iter != particleIds.end(); ++iter )
      {
      vtkIdType idx = *iter;
      particles->GetPoint(idx, pnt);
      if( mySphere->EvaluateFunction(pnt) > 0.0 )
        {
        particlesToErase.insert( idx );
        }
      } // END for all particle ids

    std::set< vtkIdType >::iterator setIter = particlesToErase.begin();
    for( ; setIter != particlesToErase.end(); ++setIter )
      {
      particleIds.remove( *setIter );
      } // END for all particle ids that are to be erased
    }

  mySphere->Delete();
}

//------------------------------------------------------------------------------
void vtkCosmoDensityProfile::ComputeDensity(vtkUnstructuredGrid* particles)
{
  int biggestSphereId = this->NumberOfSpheres-1;
  double pnt[3];
  std::list< vtkIdType > particleIds;

  for( int sphereIdx=biggestSphereId; sphereIdx >= 0; --sphereIdx )
    {
    this->FindParticlesInSphere(sphereIdx,particles,particleIds );
    this->NumParticlesInSphere[ sphereIdx ] =
        static_cast<int>( particleIds.size() );

    } // END for all spheres from bigger to smallest

  this->Controller->Barrier();

  // Reduce the sum for each sphere to process 0 which is responsible for
  // generating the final plot.
  if( this->Controller->GetNumberOfProcesses() > 1)
    {
    int *rcvBuffer = new int[this->NumberOfSpheres];
    assert("rcvBuffer is NULL!" && (rcvBuffer != NULL) );

    vtkVectorSumOperator vectorSumOperator;
    int rootProcess = 0;
    this->Controller->Reduce(
      &this->NumParticlesInSphere[0],rcvBuffer,this->NumberOfSpheres,
      &vectorSumOperator,rootProcess);

    if( this->Controller->GetLocalProcessId() == rootProcess )
      {
      for( int i=0; i < this->NumberOfSpheres; ++i )
        {
        this->NumParticlesInSphere[i] = rcvBuffer[i];
        } // END for all spheres
      } // END if rootProcess

    delete [] rcvBuffer;
    } // End if more than one process


}

//------------------------------------------------------------------------------
int vtkCosmoDensityProfile::RequestData(
    vtkInformation* vtkNotUsed(request),vtkInformationVector** inputVector,
    vtkInformationVector *outputVector )
{
  // STEP 0: Get input object
  vtkInformation *input = inputVector[0]->GetInformationObject( 0 );
  assert("pre: input information object is NULL" && (input != NULL) );
  vtkUnstructuredGrid *particles =
      vtkUnstructuredGrid::SafeDownCast(
          input->Get( vtkDataObject::DATA_OBJECT() ) );
  assert("pre: input particles is NULL!" && (particles != NULL) );

  // STEP 1: Get the 1st output
  vtkInformation *output = outputVector->GetInformationObject( 0 );
  assert( "pre: output information object is NULL" && (output != NULL) );
  vtkMultiBlockDataSet *mbds =
      vtkMultiBlockDataSet::SafeDownCast(
          output->Get( vtkDataObject::DATA_OBJECT() ) );

  // STEP 2: Get the 2nd output
  vtkInformation *output2 = outputVector->GetInformationObject( 1 );
  assert( "pre: output information object is NULL" && (output2 != NULL) );
  vtkTable *plot = vtkTable::SafeDownCast(
      output2->Get(vtkDataObject::DATA_OBJECT() ) );

  // STEP 3: Generate con-centric spheres
  this->GenerateSpheres(mbds);
  this->Controller->Barrier();

  // STEP 4: Compute density (i.e., nunmber of particles in sphere)
  this->ComputeDensity(particles);
  this->Controller->Barrier();

  if( this->Controller->GetLocalProcessId() == 0 )
    {
    vtkIdTypeArray *sphereIds       = vtkIdTypeArray::New();
    sphereIds->SetName( "Sphere ID" );
    vtkIdTypeArray *sphereDensities = vtkIdTypeArray::New();
    sphereDensities->SetName("Density");

    for( int i=0; i < this->NumberOfSpheres; ++i )
      {
      sphereIds->InsertNextValue( i );
      sphereDensities->InsertNextValue( this->NumParticlesInSphere[i] );
      } // END for all spheres

    plot->AddColumn( sphereIds );
    sphereIds->Delete();

    plot->AddColumn( sphereDensities );
    sphereDensities->Delete();
    }
  this->Controller->Barrier();

  return 1;
}
