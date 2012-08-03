/**
 * @brief An object to hold a reference to the simulation arrays for the
 * particle data, consisting of the velocity and position vectors and the
 * corresponding particle global IDs. SimulationParticles also holds a
 * reference to the halo/sub-halo tags which are computed internally.
 */
#ifndef SIMULATIONPARTICLES_H_
#define SIMULATIONPARTICLES_H_

#include "CosmologyToolsMacros.h"


namespace cosmologytools
{

class SimulationParticles
{
public:

  /**
   * @brief Default constructor
   */
  SimulationParticles();

  /**
   * @brief Destructor
   */
  ~SimulationParticles();

  /**
   * @brief Allocates the halo and subhalo arrays
   */
  void AllocateHaloAndSubHaloArrays();

  INTEGER TimeStep;     /** the current timestep */
  REAL RedShift;        /** the red-shift */
  REAL *X;              /** x-component of the particles position vector  */
  REAL *Y;              /** y-component of the particles position vector  */
  REAL *Z;              /** z-component of the particles position vector  */
  REAL *VX;             /** vx-component of the particles position vector */
  REAL *VY;             /** vy-component of the particles position vector */
  REAL *VZ;             /** vz-component of the particles position vector */
  REAL *Mass;           /** particle masses */
  REAL *Potential;      /** particle potential */
  INTEGER *GlobalIds;   /** list of global ids */

  // TODO: Are both mask & state arrays needed?
  INTEGER *Mask;        /** particle mask array */
  INTEGER *State;       /** particle state array? */

  INTEGER NumParticles; /** the total number of particles */

  INTEGER *HaloTags;    /** List of halo ids for each particle (computed) */
  INTEGER *SubHaloTags; /** List of subhalo ids for each particles (computed) */

  bool AreHalosComputed; /** Flag that indicates if the halos are computed */
  bool AreSubHalosComputed; /** Flag that indicates if subhalos are computed */

private:
  DISABLE_COPY_AND_ASSIGNMENT(SimulationParticles);
};

} /* namespace cosmogolytools */
#endif /* SIMULATIONPARTICLES_H_ */
