#include "vtkMinkowskiFilter.h"

#include <cmath>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkCell.h>
#include <vtkPolygon.h>
#include <vtkPolyhedron.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkInformation.h>
#include <vtkInformationVector.h>
#include <vtkUnstructuredGrid.h>
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
}

vtkMinkowskiFilter::~vtkMinkowskiFilter()
{
}

void vtkMinkowskiFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);
  os << indent << "vtkMinkowskiFilter\n";
}

int vtkMinkowskiFilter::RequestData(vtkInformation *vtkNotUsed(request), 
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
  // Get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
 
  // Get the input and ouptut
  vtkUnstructuredGrid *input = vtkUnstructuredGrid::SafeDownCast(
      inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::SafeDownCast(
      outInfo->Get(vtkDataObject::DATA_OBJECT()));

  output->CopyStructure(input);

  VTK_CREATE(vtkDoubleArray, S);
  VTK_CREATE(vtkDoubleArray, V);
  VTK_CREATE(vtkDoubleArray, C);
  VTK_CREATE(vtkDoubleArray, X);

  compute_mf(input, S, V, C, X);

  output->CopyStructure(input);
  output->GetPointData()->PassData(input->GetPointData());
  output->GetCellData()->PassData(input->GetCellData());
  output->GetCellData()->AddArray(S);
  output->GetCellData()->AddArray(V);
  output->GetCellData()->AddArray(C);
  output->GetCellData()->AddArray(X);

  return 1;
}

int vtkMinkowskiFilter::FillOutputPortInformation(int port, vtkInformation* info )
{
  if ( port == 0 )
  {
    info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkUnstructuredGrid" );
 
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

  vtkCellData *cell_data = ugrid->GetCellData();
  vtkFloatArray *area_array = vtkFloatArray::SafeDownCast(
      cell_data->GetArray("Areas"));
  vtkFloatArray *vol_array = vtkFloatArray::SafeDownCast(
      cell_data->GetArray("Volumes"));
  float *farea = area_array->GetPointer(0);
  float *fvol = vol_array->GetPointer(0);

  for (i=0; i<num_cells; i++)
  {
    vtkPolyhedron *cell = vtkPolyhedron::SafeDownCast(ugrid->GetCell(i));

    S_array[i] = compute_S(cell);
    V_array[i] = compute_V(ugrid, i);
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
  V->SetVoidArray(V_array, num_cells, 0);

  C->SetName("C");
  C->SetNumberOfComponents(1);
  C->SetNumberOfTuples(num_cells);
  C->SetVoidArray(C_array, num_cells, 0);

  X->SetName("X");
  X->SetNumberOfComponents(1);
  X->SetNumberOfTuples(num_cells);
  X->SetVoidArray(X_array, num_cells, 0);
  
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
  // for manually compute polyhedron volume
  return 0; // doing nothing now
}

double vtkMinkowskiFilter::compute_V(vtkUnstructuredGrid *ugrid, int cid)
{
  // for precomputed cell volume
  return ugrid->GetCellData()->GetArray("Volumes")->GetTuple1(cid);
}

double vtkMinkowskiFilter::compute_C(vtkPolyhedron *cell)
{
  return 0; // doing nothing now
}

double vtkMinkowskiFilter::compute_X(vtkPolyhedron *cell)
{
  int nfaces = cell->GetNumberOfFaces();
  int nedges = cell->GetNumberOfEdges();
  int nverts = cell->GetNumberOfPoints();

  double X = nfaces - nedges + nverts;
  
  return X;
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

  ax = fabs(N[0]);
  ay = fabs(N[1]);
  az = fabs(N[2]);

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

  return fabs(area); //area is alway positive
}

double vtkMinkowskiFilter::compute_edge_length(double v1[3], double v2[3])
{
  return sqrt(vtkMath::Dot(v1, v2));
}

double vtkMinkowskiFilter::compute_face_angle(vtkCell *f1, vtkCell *f2)
{
  return 0; // doing nothing now
}

int vtkMinkowskiFilter::get_num_edges(vtkPolyhedron *cell)
{
  return cell->GetNumberOfEdges();
}

int vtkMinkowskiFilter::get_num_faces(vtkPolyhedron *cell)
{
  return cell->GetNumberOfFaces();
}

int vtkMinkowskiFilter::get_num_verts(vtkPolyhedron *cell)
{
  return cell->GetNumberOfPoints();
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



