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


/**
 * @class MemoryLayout
 * @brief An object to hold an enum of all available memory layouts for
 * multi-dimensional arrays. Currently,two layouts are supported: (1) ROWMAJOR,
 * i.e., a c-style layout, and (2) COLMAJOR, i.e., a FORTRAN-style layout.
 */
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

// Forward declarations
class SimulationParticles;
class ForwardHaloTracker;

class CosmologyToolsManager
{
public:
  CosmologyToolsManager();
  virtual ~CosmologyToolsManager();

  // Inline class macros
  GetNSetMacro(Layout,int);
  GetNSetMacro(HaloFinder,int);
  GetNSetMacro(Communicator,MPI_Comm);
  GetNSetMacro(EnableVis,bool);
  GetNSetMacro(HaloTrackingFrequency,int);
  GetMacro(Particles,SimulationParticles*);

  /**
   * @brief Sets the particles at the given timeste/redshift.
   * @param tstep the current discrete time-step
   * @param redshift the redshift at the given time-step
   * @param px x-component of the particles position vector
   * @param py y-component of the particles position vector
   * @param pz z-component of the particles position vector
   * @param vx x-component of the particle velocity vector
   * @param vy y-component of the particle velocity vector
   * @param vz z-component of the particle velocity vector
   * @param GlobalParticlesIds the global IDs of each particle
   * @param NumberOfParticles the total number of particles
   */
  void SetParticles(
      INTEGER tstep, REAL redshift,
      REAL *px, REAL *py, REAL *pz,
      REAL *vx, REAL *vy, REAL *vz,
      INTEGER *GlobalParticlesIds,
      INTEGER NumberOfParticles);

  /**
   * @brief Uses the forward halo-tracker to track the halos at the prescribed
   * tracker frequency.
   */
  void TrackHalos();

  /**
   * @brief Calls the prescribed halo-finder to find the halos in the particle
   * dataset of the given time-step.
   */
  void FindHalos();

  /**
   * @brief Barrier synchronization with all processes.
   */
  void Barrier() {MPI_Barrier(this->Communicator);};

protected:
  int HaloFinder;
  int HaloTrackingFrequency;
  int Layout;
  MPI_Comm Communicator;
  bool EnableVis;
  SimulationParticles *Particles;

  ForwardHaloTracker *HaloTracker;

private:
  DISABLE_COPY_AND_ASSIGNMENT(CosmologyToolsManager);
};

} /* namespace cosmologytools */
#endif /* COSMOLOGYTOOLSMANAGER_H_ */
