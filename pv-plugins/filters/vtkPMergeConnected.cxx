#include "vtkPMergeConnected.h"

#include <iostream>
#include <cstring>
#include <map>

#include <vtkIdList.h>
#include <vtkPolyhedron.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
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


vtkStandardNewMacro(vtkPMergeConnected); 

vtkPMergeConnected::vtkPMergeConnected()
{
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

  this->Controller = NULL;
  this->SetController(vtkMultiProcessController::GetGlobalController());
}

vtkPMergeConnected::~vtkPMergeConnected()
{
}

void vtkPMergeConnected::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "vtkPMergeConnected filter" << "\n";
}

void vtkPMergeConnected::SetController(vtkMultiProcessController *c)
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

int vtkPMergeConnected::RequestData(vtkInformation *vtkNotUsed(request), 
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

  output->CopyStructure(input);

  int i, j, k, tc, tb;
  tb = input->GetNumberOfBlocks();
  for (i=piece; i<tb; i+=numPieces)
  {
    vtkUnstructuredGrid *ugrid = vtkUnstructuredGrid::SafeDownCast(
        input->GetBlock(i));

    VTK_CREATE(vtkUnstructuredGrid, ugrid_out);
    ugrid_out->SetPoints(ugrid->GetPoints());

    vtkCellData   *cd = ugrid->GetCellData();
    vtkPointData  *pd = ugrid->GetPointData();
    vtkCellData  *ocd = ugrid_out->GetCellData();
    vtkPointData *opd = ugrid_out->GetPointData();

    ocd->CopyStructure(cd);
    opd->CopyStructure(pd);

    vtkIdTypeArray *prid_array = vtkIdTypeArray::SafeDownCast(
        pd->GetArray("RegionId"));
    vtkIdTypeArray *crid_array = vtkIdTypeArray::SafeDownCast(
        cd->GetArray("RegionId"));
    vtkIdTypeArray *oprid_array = vtkIdTypeArray::SafeDownCast(
        opd->GetArray("RegionId"));
    vtkIdTypeArray *ocrid_array = vtkIdTypeArray::SafeDownCast(
        ocd->GetArray("RegionId"));

    //oprid_array->SetNumberOfComponents(1);
    //ocrid_array->SetNumberOfComponents(1);

    //checking RegionId range
    double rid_range[2];
    crid_array->GetRange(rid_range, 0);
    int num_regions = rid_range[1] - rid_range[0] + 1;

    //initialize the size of cell data arrays
    for (j=0; j<cd->GetNumberOfArrays(); j++)
      ocd->GetArray(j)->SetNumberOfTuples(num_regions);

    tc = ugrid->GetNumberOfCells();
    int new_id = 0;
    for (j=rid_range[0]; j<=rid_range[1]; j++)
    {
      //compute face stream of merged polyhedron cell
      vtkIdList *mcell = MergeCellsOnRegionId(ugrid, j);
      ugrid_out->InsertNextCell(VTK_POLYHEDRON, mcell);
    
      //"RegionId" cells keep their old id
      ocrid_array->SetValue(new_id, j);
      for (k=0; k<cd->GetNumberOfArrays(); k++)
      {
        vtkDataArray *array = cd->GetArray(k);
        vtkDataArray *oarray = ocd->GetArray(k);
        if (!std::strcmp(array->GetName(), "Volume"))
        {
          vtkFloatArray *vol_array = vtkFloatArray::SafeDownCast(array);
          float vol = MergeCellDataOnRegionId(vol_array, crid_array, j);
          vtkFloatArray *ovol_array = vtkFloatArray::SafeDownCast(oarray);
          ovol_array->SetValue(new_id, vol);
        }
      }
      new_id ++;
    }

    output->SetBlock(i, ugrid_out);
  }

  return 1;
}

std::string vtkPMergeConnected::IdsToKeystr(vtkIdList* ids)
{
  return "key";
}

struct cmp_ids
{
  bool operator()(char const *a, char const *b)
  {
    int i, ret = 0;
    vtkIdList *list_a = (vtkIdList *)a;
    vtkIdList *list_b = (vtkIdList *)b;
    
    if (list_a->GetNumberOfIds() == list_b->GetNumberOfIds())
    {
      for (i=0; i<list_a->GetNumberOfIds(); i++)
      {
        if (list_a->GetId(i) != list_b->GetId(i))
        {
          ret = list_a->GetId(i) - list_b->GetId(i);
          break;
        } 
      }
    }
    else
      ret = list_a->GetNumberOfIds() - list_b->GetNumberOfIds();

    return ret;
  }
};


//face stream of a polyhedron cell in the following format: 
//numCellFaces, numFace0Pts, id1, id2, id3, numFace1Pts,id1, id2, id3, ...
vtkIdList* vtkPMergeConnected::MergeCellsOnRegionId(vtkUnstructuredGrid *ugrid, int target)
{
  int i, j;

  VTK_CREATE(vtkIdList, facestream);
  vtkIdTypeArray *rids = vtkIdTypeArray::SafeDownCast(ugrid->GetCellData()->GetArray("RegionId"));

  //initially set -1 for number of faces in facestream
  facestream->InsertNextId(-1);

  std::map<vtkIdList*, int> face_map;

  for (i=0; i<rids->GetNumberOfTuples(); i++)
  {
    if (rids->GetTuple1(i) == target)
    {
      vtkPolyhedron *cell = vtkPolyhedron::SafeDownCast(ugrid->GetCell(i));
      for (j=0; j<cell->GetNumberOfFaces(); j++)
      {
        vtkCell *face = cell->GetFace(j);
        vtkIdList *pts = face->GetPointIds();
        std::string key = IdsToKeyStr(pts);
        
        std::map<std::string, int>::iterator it = face_map.find(key);
        if (it == edge_map.end())
          it->second = 0;
        else
          it->second ++;
      }
    }
  }
}

float vtkPMergeConnected::MergeCellDataOnRegionId(vtkFloatArray *data_array, vtkIdTypeArray *rid_array, vtkIdType target)
{
  int i;
  float val = 0;

  for (i=0; i<rid_array->GetNumberOfTuples(); i++)
    if (rid_array->GetTuple1(i) == target)
      val += data_array->GetTuple1(i);

  return val;
}

int vtkPMergeConnected::FillOutputPortInformation(int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkMultiBlockDataSet" );
 
    return 1;
  }
 
  return 0;
}

