#include "VoronoiFilter.h"

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkFieldData.h>
#include <vtkFloatArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiProcessController.h>
#include <vtkObjectFactory.h>
#include <vtkSmartPointer.h>
#include <vtkStreamingDemandDrivenPipeline.h>
#include <vtkUnstructuredGrid.h>

#define VTK_CREATE(type, name) \
  vtkSmartPointer<type> name = vtkSmartPointer<type>::New()
#define VTK_NEW(type, name) \
  name = vtkSmartPointer<type>::New()


vtkStandardNewMacro(VoronoiFilter);

VoronoiFilter::VoronoiFilter()
{
  this->VolumeRange[0] = 0.0;
  this->VolumeRange[0] = 1.0;

  this->AreaRange[0] = 0.0;
  this->AreaRange[0] = 1.0;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

VoronoiFilter::~VoronoiFilter()
{
}

void VoronoiFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "Volume Range: ("
     << this->VolumeRange[0] << ", "
     << this->VolumeRange[1] << ")\n";
  os << indent << "Area Range: ("
     << this->AreaRange[0] << ", "
     << this->AreaRange[1] << ")\n";
}

void VoronoiFilter::SetController(vtkMultiProcessController *c)
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

int VoronoiFilter::RequestData(vtkInformation *vtkNotUsed(request),
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

  if (sum > 0) //from example, not sure whether to happen
  {
    return 1; //to be changed to the right error code
  }

  if (!contr) //from example, not sure whether to happen
  {
    return 1;
  }

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

    //output
    VTK_CREATE(vtkUnstructuredGrid, ugrid_out);
    VTK_CREATE(vtkFloatArray, area_array_out);
    area_array_out->SetName("Areas");
    area_array_out->SetNumberOfComponents(1);
    VTK_CREATE(vtkFloatArray, vol_array_out);
    vol_array_out->SetName("Volumes");
    vol_array_out->SetNumberOfComponents(1);

    ugrid_out->SetPoints(ugrid->GetPoints());
    for (j=0; j<tc; j++)
    {
      if (area_array->GetValue(j) >= AreaRange[0] &&
          area_array->GetValue(j) <= AreaRange[1] &&
          vol_array->GetValue(j) >= VolumeRange[0] &&
          vol_array->GetValue(j) <= VolumeRange[1])
      {
        VTK_CREATE(vtkIdList, ptIds);
        ugrid->GetFaceStream(j, ptIds);
        ugrid_out->InsertNextCell(VTK_POLYHEDRON, ptIds);

        area_array_out->InsertNextTuple1(area_array->GetValue(j));
        vol_array_out->InsertNextTuple1(vol_array->GetValue(j));
      }
    }
    ugrid_out->GetCellData()->AddArray(area_array_out);
    ugrid_out->GetCellData()->AddArray(vol_array_out);

    output->SetBlock(i, ugrid_out);
  }

  return 1;
}

int VoronoiFilter::FillOutputPortInformation(int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );

    return 1;
  }

  return 0;
}
