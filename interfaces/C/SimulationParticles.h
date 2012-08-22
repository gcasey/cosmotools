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

  INTEGER TimeStep;         /** the current timestep */
  REAL RedShift;            /** the red-shift */
  POSVEL_T *X;              /** x-component of the particles position vector  */
  POSVEL_T *Y;              /** y-component of the particles position vector  */
  POSVEL_T *Z;              /** z-component of the particles position vector  */
  POSVEL_T *VX;             /** vx-component of the particles position vector */
  POSVEL_T *VY;             /** vy-component of the particles position vector */
  POSVEL_T *VZ;             /** vz-component of the particles position vector */
  REAL *Mass;               /** particle masses */
  POTENTIAL_T *Potential;   /** particle potential */
  ID_T *GlobalIds;          /** list of global ids */

  // TODO: Are both mask & state arrays needed?
  MASK_T *Mask;        /** particle mask array */
  STATUS_T *State;       /** particle state array? */

  INTEGER NumParticles; /** the total number of particles */

  INTEGER *HaloTags;    /** List of halo ids for each particle (computed) */
  INTEGER *SubHaloTags; /** List of subhalo ids for each particles (computed)*/

  bool AreHalosComputed; /** Flag that indicates if the halos are computed */
  bool AreSubHalosComputed; /** Flag that indicates if sub-halos are computed*/

private:
  DISABLE_COPY_AND_ASSIGNMENT(SimulationParticles);
};

} /* namespace cosmogolytools */
#endif /* SIMULATIONPARTICLES_H_ */
