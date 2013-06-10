/**
 * @brief An object to hold a reference to the simulation arrays for the
 * particle data, consisting of the velocity and position vectors and the
 * corresponding particle global IDs. SimulationParticles also holds a
 * reference to the halo/sub-halo tags which are computed internally.
 */
#ifndef SIMULATIONPARTICLES_H_
#define SIMULATIONPARTICLES_H_

#include "CosmoToolsMacros.h"


namespace cosmotk
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
    * @brief Sets the particle information at the given time-step.
    * @param tstep the current time-step
    * @param redshift the redshift at the given time-step
    * @param px the x-coordinate of the particle position vector
    * @param py the y-coordinate of the particle position vector
    * @param pz the z-coordinate of the particle position vector
    * @param vx the x-coordinate of the particle velocity vector
    * @param vy the y-coordinate of the particle velocity vector
    * @param vz the z-coordinate of the particle velocity vector
    * @param mass array of particle masses
    * @param potential array of particle potential
    * @param tags array of particle tags
    * @param mask array of particle masking
    * @param status array of particle status
    * @param N the total number of particles for the given process
    */
  void SetParticles(
      INTEGER tstep, REAL redShift,
      POSVEL_T *px, POSVEL_T *py, POSVEL_T *pz,
      POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz,
      POSVEL_T *mass, POTENTIAL_T *potential,
      ID_T *tags, MASK_T *mask, STATUS_T *status,
      ID_T N);

  INTEGER TimeStep;       /** the current timestep */
  REAL RedShift;          /** the red-shift */
  POSVEL_T *X;            /** x-component of the particles position vector  */
  POSVEL_T *Y;            /** y-component of the particles position vector  */
  POSVEL_T *Z;            /** z-component of the particles position vector  */
  POSVEL_T *VX;           /** vx-component of the particles position vector */
  POSVEL_T *VY;           /** vy-component of the particles position vector */
  POSVEL_T *VZ;           /** vz-component of the particles position vector */
  POSVEL_T *Mass;         /** particle masses */
  POTENTIAL_T *Potential; /** particle potential */
  ID_T *GlobalIds;        /** list of global ids */

  // TODO: Are both mask & state arrays needed?
  MASK_T   *Mask;        /** particle mask array */
  STATUS_T *State;       /** particle state array? */

  ID_T NumParticles; /** the total number of particles */

private:
  DISABLE_COPY_AND_ASSIGNMENT(SimulationParticles);
};

} // END cosmotk namespace
#endif /* SIMULATIONPARTICLES_H_ */
