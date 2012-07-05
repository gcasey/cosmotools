#ifndef __vtkMinkowskiFilter_h
#define __vtkMinkowskiFilter_h

#include "vtkUnstructuredGridAlgorithm.h"

class vtkUnstructuredGrid;
class vtkCell;
class vtkDoubleArray;
class vtkPolyhedron;

class VTK_EXPORT vtkMinkowskiFilter : public vtkUnstructuredGridAlgorithm 
{
 public:
  static vtkMinkowskiFilter *New();
  vtkTypeMacro(vtkMinkowskiFilter, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

 protected:
  vtkMinkowskiFilter();
  ~vtkMinkowskiFilter();

  int RequestData(vtkInformation*,
                  vtkInformationVector**,
                  vtkInformationVector*);
  int FillOutputPortInformation(int, vtkInformation*);

 private:
  vtkMinkowskiFilter(const vtkMinkowskiFilter&);  // Not implemented.
  void operator=(const vtkMinkowskiFilter&);  // Not implemented.

  void compute_mf(vtkUnstructuredGrid *ugrid, vtkDoubleArray *S, 
      vtkDoubleArray *V, vtkDoubleArray *C, vtkDoubleArray *X);
  double compute_S(vtkPolyhedron *cell); //surface area
  double compute_V(vtkPolyhedron *cell); //volume1
  double compute_V(vtkUnstructuredGrid *ugrid, int cid); //volume2
  double compute_C(vtkPolyhedron *cell); //integrated mean curvature
  double compute_X(vtkPolyhedron *cell); //euler characteristic

  void compute_normal(vtkCell *face, double normal[3]);
  double compute_face_area(vtkCell *face);
  double compute_edge_length(double v1[3], double v2[3]);
  double compute_face_angle(vtkCell *f1, vtkCell *f2);
  int get_num_edges(vtkPolyhedron *cell);
  int get_num_faces(vtkPolyhedron *cell);
  int get_num_verts(vtkPolyhedron *cell);
};

#endif //  __vtkMinkowskiFilter_h

