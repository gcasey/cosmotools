/**
 * @brief C interface for cosmology tools
 */
#ifndef COSMOLOGY_TOOLS_H_
#define COSMOLOGY_TOOLS_H_

#include "CosmologyToolsAPIMangling.h" // auto-generated
#include "CosmologyToolsManager.h"
#include <mpi.h>

using namespace cosmologytools;

CosmologyToolsManager *CosmoToolsManager;


#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initialize the cosmology tools interface with an MPI communicator.
 * @param comm a C MPI communicator that will be used.
 */
void CosmologyInit(MPI_Comm *comm);

/**
 * @brief Initialize the cosmology tool with the given Fortran MPI communicator.
 * @param fcomm a Fortran MPI communicator that will be used.
 * @note The fortran communicator is converted to the corresponding C
 * communicator internally.
 */
void CosmologyFinit(MPI_Fint *fcomm);

/**
 * @brief Enable/Disable in-situ visualization.
 * @pre CosmoToolsManager != NULL.
 */
void CosmologyEnableVis();
void CosmologyDisableVis();

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
void CosmologySetParticles(
      INTEGER *tstep, REAL *redshift,
      REAL *px, REAL *py, REAL *pz,
      REAL *vx, REAL *vy, REAL *vz,
      INTEGER *GlobalParticlesIds,
      INTEGER *NumberOfParticles);

/**
 * @brief Sets the memory layout that the interface will be using.
 * @param memoryLayout the memory layout
 * @note The memory layout is defined as follows:
 * <ul>
 *   <li>ROWMAJOR, i.e., xyz xyz xyz...</li>
 *   <li>COLMAJOR, i.e., xxx...yyy...zzz</li>
 * </ul>
 */
void CosmologySetMemoryLayout(int* memoryLayout);

/**
 * @brief Sets the analysis time-steps the cosmology tools executes
 * @param tsteps user-supplied array of time-steps
 * @param N the number of time-steps, i.e., size of the array.
 */
void CosmologySetAnalysisTimeSteps(int *tsteps, int N);

/**
 * @brief Sets the frequency at which the tracker will be invoked.
 * @param frequency the frequency for the halo tracker
 */
void CosmologySetTrackerFrequency(int *frequency);

/**
 * @brief Fills the user-supplied buffer with the halo tags for each particle.
 * @param haloTags user-supplied buffer to store the halo tags
 * @note particles that are not inside a halo will have a value of -1.
 */
void CosmologyGetHaloIds(INTEGER *haloTags);

/**
 * @brief Fills the user-supplied buffer with the subhalo tags for each particle.
 * @param subHaloTags user-supplied buffer to store the sub-halo tags
 * @note particles that aere not inside a sub-halo will have a valu of -1.
 */
void CosmologyGetSubHaloIds(INTEGER *subHaloTags);

/**
 * @brief Sets the halo-finder to use
 * @param haloFinder the halo-finder to use
 * @see HaloFinders.h for a list of halofinders to use.
 */
void CosmologySetHaloFinder(int *haloFinder);

/**
 * @brief Uses the forward halo-tracker to track halos at the prescribed
 * tracker frequency.
 * @note This method is used in combination with CosmologySetParticles
 */
void CosmologyTrackHalos();

/**
 * @brief Calls the prescribed halo finder to find the halos in the particle
 * dataset of the given time-step.
 * @note This method is used in combination with CosmologySetParticles
 */
void CosmologyFindHalos();

/**
 * @brief Finalize the cosmology tools interface.
 * @note Blocks until all processes call this method.
 * @pre CosmoToolsManager != NULL
 */
void CosmologyFinalize();

#ifdef __cplusplus
}
#endif

#endif /* COSMOLOGY_TOOLS_H_ */
