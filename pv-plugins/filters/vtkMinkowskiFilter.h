#ifndef __vtkMinkowskiFilter_h
#define __vtkMinkowskiFilter_h

#include "vtkMultiBlockDataSetAlgorithm.h"

class vtkMultiProcessController;
class vtkUnstructuredGrid;

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

  void compute_mf(vtkUnstructuredGrid *ugrid, 
                  double &S, 
                  double &V, 
                  double &C, 
                  double &X);
  double compute_S(vtkUnstructuredGrid *ugrid); //surface area
  double compute_V(vtkUnstructuredGrid *ugrid); //volume
  double compute_C(vtkUnstructuredGrid *ugrid); //integrated mean curvature
  double compute_X(vtkUnstructuredGrid *ugrid); //euler characteristic

  double compute_face_area(vtkCell *face, vtkPoints *pts);
  double compute_edge_length(vtkPoint);
  double compute_face_angle(vtkCell *f1, vtkCell *f2);
  double compute_get_num_edges(vtkCell *cell);
};

#endif //  __vtkMinkowskiFilter_h

