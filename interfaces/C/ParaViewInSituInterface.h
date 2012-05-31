/**
 * @brief A simple interface that the simulation/co-processor can use for
 * interacting with ParaView.
 */
#ifndef PARAVIEWINSITUINTERFACE_H_
#define PARAVIEWINSITUINTERFACE_H_

#include "CosmologyToolsMacros.h"

// Forward declarations
class vtkLiveInsituLink;

namespace cosmologytools
{

// Forward declarations within cosmologytools namespace
class SimulationParticles;

class ParaViewInSituInterface
{
public:
  ParaViewInSituInterface();
  virtual ~ParaViewInSituInterface();

  // In-line methods
  GetMacro(Port,int);

  /**
   * @brief Initialize the in-situ paraview interface
   * @param port optionally provide a port number, by default 2222 is used.
   */
  void Initialize(int port=2222);

  /**
   * @brief Update particles for the current red-shift.
   * @param particles pointer to the SimulationParticles data-structure.
   * @pre particles != NULL
   */
  void Update(SimulationParticles *particles);

  /**
   * @brief Finalize interface
   */
  void Finalize();

protected:
  int Port;
  bool IsInitialized;
  vtkLiveInsituLink *PVLink;

private:
  DISABLE_COPY_AND_ASSIGNMENT(ParaViewInSituInterface);
};

} /* namespace cosmologytools */
#endif /* PARAVIEWINSITUINTERFACE_H_ */
