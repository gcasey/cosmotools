#ifndef __vtkMinkowskiFilter_h
#define __vtkMinkowskiFilter_h

#include "vtkMultiBlockDataSetAlgorithm.h"

class vtkMultiProcessController;

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
};

#endif //  __vtkMinkowskiFilter_h

