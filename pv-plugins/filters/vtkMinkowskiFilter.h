#ifndef __vtkMinkowskiFilter_h
#define __vtkMinkowskiFilter_h

#include "vtkMultiBlockDataSetAlgorithm.h"

class vtkMultiProcessController;
class vtkUnstructuredGrid;
class vtkCell;
class vtkDoubleArray;
class vtkPoints;
class vtkPolyhedron;

class VTK_EXPORT vtkMinkowskiFilter : public vtkMultiBlockDataSetAlgorithm 
{
 public:
  static vtkMinkowskiFilter *New();
  vtkTypeMacro(vtkMinkowskiFilter, vtkMultiBlockDataSetAlgorithm);
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

  int NumProcesses;
  int MyId;
  vtkMultiProcessController *Controller;
  void SetController(vtkMultiProcessController *c);

  void compute_mf(vtkUnstructuredGrid *ugrid, vtkDoubleArray *S, 
      vtkDoubleArray *V, vtkDoubleArray *C, vtkDoubleArray *X);
  double compute_S(vtkPolyhedron *cell); //surface area
  double compute_V(vtkPolyhedron *cell); //volume
  double compute_C(vtkPolyhedron *cell); //integrated mean curvature
  double compute_X(vtkPolyhedron *cell); //euler characteristic

  void compute_normal(vtkCell *face, double normal[3]);
  double compute_face_area(vtkCell *face);
  double compute_edge_length(double v1[3], double v2[3]);
  double compute_face_angle(vtkCell *f1, vtkCell *f2);
  int get_num_edges(vtkCell *cell);
};

#endif //  __vtkMinkowskiFilter_h

