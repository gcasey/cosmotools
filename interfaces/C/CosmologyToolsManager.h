/**
 * @brief CosmologyToolsManager to keep persistent information through the
 * life-time of the simulation.
 */
#ifndef COSMOLOGYTOOLSMANAGER_H_
#define COSMOLOGYTOOLSMANAGER_H_

#include "CosmologyToolsMacros.h"
#include <mpi.h>

namespace cosmologytools
{

class MemoryLayout
{
public:
  enum
    {
    ROWMAJOR=0,
    COLMAJOR=1,
    NUMBER_OF_LAYOUTS
    };
};

class CosmologyToolsManager
{
public:
  CosmologyToolsManager();
  virtual ~CosmologyToolsManager();

  // Inline class macros
  GetNSetMacro(Layout,int);
  GetNSetMacro(Communicator,MPI_Comm);
  GetNSetMacro(EnableVis,bool);
  GetNSetMacro(ArrayStride,int);
  GetNSetMacro(HaloTrackingFrequency,int);

  /**
   * @brief Barrier synchronization with all processes.
   */
  void Barrier() {MPI_Barrier(this->Communicator);};

protected:
  int HaloFinder;
  int HaloTrackingFrequency;
  int Layout;
  int ArrayStride;
  MPI_Comm Communicator;
  bool EnableVis;

private:
  DISABLE_COPY_AND_ASSIGNMENT(CosmologyToolsManager);
};

} /* namespace cosmologytools */
#endif /* COSMOLOGYTOOLSMANAGER_H_ */
