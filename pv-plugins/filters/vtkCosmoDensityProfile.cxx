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

#include <cassert>

vtkStandardNewMacro(vtkCosmoDensityProfile);

//------------------------------------------------------------------------------
vtkCosmoDensityProfile::vtkCosmoDensityProfile()
{
  this->Controller = vtkMultiProcessController::GetGlobalController();
  this->Center[0]  = this->Center[1] = this->Center[2] = 0.0;
  this->Radius     = 5.0;
  this->SetNumberOfInputPorts( 1 );
  this->SetNumberOfOutputPorts( 2 );
}

//------------------------------------------------------------------------------
vtkCosmoDensityProfile::~vtkCosmoDensityProfile()
{
  // TODO Auto-generated destructor stub
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
void vtkCosmoDensityProfile::GenerateSpheres()
{

}

//------------------------------------------------------------------------------
int vtkCosmoDensityProfile::RequestData(
    vtkInformation* vtkNotUsed(request),vtkInformationVector** inputVector,
    vtkInformationVector *outputVector )
{
  // STEP 0: Get input object
  vtkInformation *input = inputVector[0]->GetInformationObject( 0 );
  assert("pre: input information object is NULL" && (input != NULL) );
//  vtkUnstructuredGrid *particles =
//      vtkUnstructuredGrid::SafeDowncast(
//          input->Get( vtkDataObject::DATA_OBJECT() ) );
//  assert("pre: input particles is NULL!" && (particles != NULL) );

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

  // TODO: implement this
  this->Controller->Barrier();
  return 1;
}
