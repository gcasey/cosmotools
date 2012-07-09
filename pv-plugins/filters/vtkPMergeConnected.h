#ifndef __vtkPMergeConnected_h
#define __vtkPMergeConnected_h

#include "vtkMultiBlockDataSetAlgorithm.h"

class vtkMultiProcessController;
class vtkUnstructuredGrid;
class vtkIdList;
class vtkFloatArray;
class vtkIdTypeArray;

class VTK_EXPORT vtkPMergeConnected : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkPMergeConnected* New();
  vtkTypeMacro(vtkPMergeConnected, vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  struct FaceWithKey
  {
    int num_pts;
    vtkIdType *key, *orig;
  };
  struct cmp_ids;

protected:
  vtkPMergeConnected();
  ~vtkPMergeConnected();

  int RequestData(vtkInformation*, vtkInformationVector**, vtkInformationVector*);
  int FillOutputPortInformation(int port, vtkInformation* info);

private:
  vtkPMergeConnected(const vtkPMergeConnected&);  // Not implemented.
  void operator=(const vtkPMergeConnected&);  // Not implemented.

  //parallelism
  int NumProcesses;
  int MyId;
  vtkMultiProcessController *Controller;
  void SetController(vtkMultiProcessController *c);

  //filter
  vtkIdList* MergeCellsOnRegionId(vtkUnstructuredGrid *ugrid, int target);
  float MergeCellDataOnRegionId(vtkFloatArray *data_array, vtkIdTypeArray *rid_array, vtkIdType target);

  void delete_key(FaceWithKey *key);
  FaceWithKey* IdsToKey(vtkIdList* ids);
};

#endif
