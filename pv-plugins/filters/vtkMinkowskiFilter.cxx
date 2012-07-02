#include "vtkMinkowskiFilter.h"

#include <cmath>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkPolyhedron.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
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

void vtkMinkowskiFilter::compute_mf(vtkUnstructuredGrid *ugrid, 
    vtkDoubleArray *S, vtkDoubleArray *V, vtkDoubleArray *C, 
    vtkDoubleArray *X)
{
  int i;
  int num_cells = ugrid->GetNumberOfCells();
  double *S_array, *V_array, *C_array, *X_array;

  S_array = new double[num_cells];
  V_array = new double[num_cells];
  C_array = new double[num_cells];
  X_array = new double[num_cells];

  vtkIdType *ptIds;
  vtkIdType nfaces;
  for (i=0; i<num_cells; i++)
  {
    ugrid->GetFaceStream(i, nfaces, ptIds);
    VTK_CREATE(vtkPolyhedron, cell);
    cell->SetFaces(ptIds);

    S_array[i] = compute_S(cell);
    V_array[i] = compute_V(cell);
    C_array[i] = compute_C(cell);
    X_array[i] = compute_X(cell);
  }

  S->SetName("S");
  S->SetNumberOfComponents(1);
  S->SetNumberOfTuples(num_cells);
  S->SetVoidArray(S_array, num_cells, 0);

  V->SetName("V");
  V->SetNumberOfComponents(1);
  V->SetNumberOfTuples(num_cells);
  V->SetVoidArray(S_array, num_cells, 0);

  C->SetName("C");
  C->SetNumberOfComponents(1);
  C->SetNumberOfTuples(num_cells);
  C->SetVoidArray(S_array, num_cells, 0);

  X->SetName("X");
  X->SetNumberOfComponents(1);
  X->SetNumberOfTuples(num_cells);
  X->SetVoidArray(S_array, num_cells, 0);
  
}

double vtkMinkowskiFilter::compute_S(vtkPolyhedron *cell)
{
  int i;
  int num_faces = cell->GetNumberOfFaces();
  double area = 0.0;
  
  vtkCell *face;
  for (i=0; i<num_faces; i++)
  {
    face = cell->GetFace(i);
    area += compute_face_area(face);
  }

  return area;
}

double vtkMinkowskiFilter::compute_V(vtkPolyhedron *cell)
{
  return 0; // doing nothing now
}

double vtkMinkowskiFilter::compute_C(vtkPolyhedron *cell)
{
  return 0; // doing nothing now
}

double vtkMinkowskiFilter::compute_X(vtkPolyhedron *cell)
{
  return 0; // doing nothing now
}

double vtkMinkowskiFilter::compute_face_area(vtkCell *face)
{
  int i, j, k;
  int coord; // coord to ignore: 1=x, 2=y, 3=z
  double area = 0;
  double an, ax, ay, az, N[3];

  int num_verts = face->GetNumberOfPoints();
  vtkPoints *verts = face->GetPoints();
  compute_normal(face, N);

  ax = abs(N[0]);
  ay = abs(N[1]);
  az = abs(N[2]);

  coord = 3;
  if (ax > ay)
  {
    if (ax > az) 
      coord = 1;
  }
  else if (ay > az)
    coord = 2;

  // compute area of the 2D projection
  double v1[3], v2[3], v3[3];
  for (i=1, j=2, k=0; i<=num_verts; i++, j++, k++)
  {
    verts->GetPoint(i % num_verts, v1);
    verts->GetPoint(j % num_verts, v2);
    verts->GetPoint(k % num_verts, v3);

    if (coord == 1)
      area += v1[1] * (v2[2] - v3[2]);
    else if (coord == 2)
      area += v1[0] * (v2[2] - v3[2]);
    else if (coord == 3)
      area += v1[0] * (v2[1] - v3[1]);
  }

  // scale to get area before projection
  an = sqrt(ax * ax + ay * ay + az * az);
  if (coord == 1)
    area *= an / (2 * ax);
  else if (coord == 2)
    area *= an / (2 * ay);
  else if (coord == 3)
    area *= an / (2 * az);

  return area;
}

double vtkMinkowskiFilter::compute_edge_length(double v1[3], double v2[3])
{
  return sqrt(vtkMath::Dot(v1, v2));
}

double vtkMinkowskiFilter::compute_face_angle(vtkCell *f1, vtkCell *f2)
{
  return 0; // doing nothing now
}

int vtkMinkowskiFilter::get_num_edges(vtkCell *cell)
{
  return 0; // doing nothing now
}

//compute normal of a face using Newell's method
void vtkMinkowskiFilter::compute_normal(vtkCell *face, double normal[3])
{
  int i;
  int num_verts = face->GetNumberOfPoints();
  vtkPoints *verts = face->GetPoints();
  double vcurr[3], vnext[3];

  normal[0] = 0.0;
  normal[1] = 0.0;
  normal[2] = 0.0;

  int curr, next;
  for (i = 0; i < num_verts; i++)
  {
    curr = i;
    next = (i + 1) % num_verts;
    verts->GetPoint(curr, vcurr);
    verts->GetPoint(next, vnext);
    normal[0] += (vcurr[1] - vnext[1]) * (vcurr[2] + vnext[2]);
    normal[1] += (vcurr[2] - vnext[2]) * (vcurr[0] + vnext[0]);
    normal[2] += (vcurr[0] - vnext[0]) * (vcurr[1] + vnext[1]);
  }

  double mag = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
		   normal[2] * normal[2]);
  // normalize
  normal[0] /= mag;
  normal[1] /= mag;
  normal[2] /= mag;

  // direction is inward, need to invert
  normal[0] *= -1.0;
  normal[1] *= -1.0;
  normal[2] *= -1.0;
}



