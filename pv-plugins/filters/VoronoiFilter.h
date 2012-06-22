#ifndef __VoronoiFilter_h
#define __VornoniFilter_h

#include "vtkMultiBlockDataSetAlgorithm.h"

class vtkMultiProcessController;

class VTK_EXPORT VoronoiFilter : public vtkMultiBlockDataSetAlgorithm
{
public:
  static VoronoiFilter* New();
  vtkTypeMacro(VoronoiFilter, vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  void SetMinVolume(double val) { VolumeRange[0] = val; }
  double GetMinVolume() { return VolumeRange[0]; }
  void SetMaxVolume(double val) { VolumeRange[1] = val; }
  double GetMaxVolume() { return VolumeRange[1]; }
  vtkSetVector2Macro(VolumeRange, double);
  vtkGetVectorMacro(VolumeRange, double, 2);

  void SetMinArea(double val) { AreaRange[0] = val; }
  double GetMinArea() { return AreaRange[0]; }
  void SetMaxArea(double val) { AreaRange[1] = val; }
  double GetMaxArea() { return AreaRange[1]; }
  vtkSetVector2Macro(AreaRange, double);
  vtkGetVectorMacro(AreaRange, double, 2);

protected:
  VoronoiFilter();
  ~VoronoiFilter();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int FillOutputPortInformation(int port, vtkInformation* info);

  double VolumeRange[2];
  double AreaRange[2];

private:
  VoronoiFilter(const VoronoiFilter&);  // Not implemented.
  void operator=(const VoronoiFilter&);  // Not implemented.

  //parallelism
  int NumProcesses;
  int MyId;
  vtkMultiProcessController *Controller;
  void SetController(vtkMultiProcessController *c);

  //filter

};

#endif
