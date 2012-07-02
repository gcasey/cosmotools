#include "vtkMinkowskiFilter.h"

#include "vtkMinkowskiFilter.h"

#include <cmath>
#include <vtkMath.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiProcessController.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStreamingDemandDrivenPipeline.h>

#include <vtkSmartPointer.h>
#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
#define VTK_NEW(type, name) \
  name = vtkSmartPointer<type>::New()

vtkStandardNewMacro(vtkMinkowskiFilter); 

vtkMinkowskiFilter::vtkMinkowskiFilter()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

vtkMinkowskiFilter::~vtkMinkowskiFilter()
{
}

void vtkMinkowskiFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "vtkMinkowskiFilter\n";
}

void vtkMinkowskiFilter::SetController(vtkMultiProcessController *c)
{
  if ((c == NULL) || (c->GetNumberOfProcesses() == 0))
  {
    this->NumProcesses = 1;
    this->MyId = 0;
  }

  if (this->Controller == c)
  {
    return;
  }

  this->Modified();

  if (this->Controller != NULL)
  {
    this->Controller->UnRegister(this);
    this->Controller = NULL;
  }

  if (c == NULL)
  {
    return;
  }

  this->Controller = c;

  c->Register(this);
  this->NumProcesses = c->GetNumberOfProcesses();
  this->MyId = c->GetLocalProcessId();
}

int vtkMinkowskiFilter::RequestData(vtkInformation *vtkNotUsed(request), 
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
 
  // Get the input and ouptut
  vtkMultiBlockDataSet *input = vtkMultiBlockDataSet::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkMultiBlockDataSet *output = vtkMultiBlockDataSet::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  int piece, numPieces;
  piece = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER());
  numPieces = outInfo->Get(
      vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES());

  //read file
  vtkMultiProcessController *contr = this->Controller;

  int sum = 0;
  int oops = ((piece != this->MyId) || (numPieces != this->NumProcesses));

  contr->Reduce(&oops, &sum, 1, vtkCommunicator::SUM_OP, 0);
  contr->Broadcast(&sum, 1, 0);

  if (sum > 0) 
    return 1; 

  if (!contr)
    return 1;

  int i, j, tc, tb;
  tb = input->GetNumberOfBlocks();
  output->SetNumberOfBlocks(tb);

  for (i=piece; i<tb; i+=numPieces)
  {
    //input
    vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::SafeDownCast(
        input->GetBlock(i));
    vtkCellData *cell_data = ugrid->GetCellData();
    vtkFloatArray *area_array = vtkFloatArray::SafeDownCast(
        cell_data->GetArray("Areas"));
    vtkFloatArray *vol_array = vtkFloatArray::SafeDownCast(
        cell_data->GetArray("Volumes"));
    tc = ugrid->GetNumberOfCells();
  }

  return 1;
}

int vtkMinkowskiFilter::FillOutputPortInformation(int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
 
    return 1;
  }
 
  return 0;
}

void vtkMinkowskiFilter::compute_mf(vtkUnstructuredGrid *ugrid, double &S, double &V, double &C, double &X)
{
  S = compute_S(ugrid);
  V = compute_V(ugrid);
  C = compute_C(ugrid);
  X = compute_X(ugrid);
}

double vtkMinkowskiFilter::compute_S(vtkUnstructuredGrid *ugrid)
{
  
}

double vtkMinkowskiFilter::compute_V(vtkUnstructuredGrid *ugrid)
{

}

double vtkMinkowskiFilter::compute_C(vtkUnstructuredGrid *ugrid)
{

}

double vtkMinkowskiFilter::compute_X(vtkUnstructuredGrid *ugrid)
{

}

double vtkMinkowskiFilter::compute_face_area(vtkCell *face)
{
  int i, j;
  double area = 0;

  j = num_points-1;  // The last vertex is the 'previous' one to the first

  for (i=0; i<num_points; i++)
    { area = area +  (X[j]+X[i]) * (Y[j]-Y[i]); 
      j = i;  //j is previous vertex to i
    }
  return area/2;
  
}

double vtkMinkowskiFilter::compute_edge_length(double *v1, double *v2)
{
  return sqrt(vtkMath::Dot(v1, v2));
}

double vtkMinkowskiFilter::compute_face_angle(vtkCell *f1, vtkCell *f2)
{

}

double vtkMinkowskiFilter::compute_get_num_edges(vtkCell *cell)
{

}


