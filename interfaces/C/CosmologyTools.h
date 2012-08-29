/**
 * @brief C interface for cosmology tools
 */
#ifndef COSMOLOGY_TOOLS_H_
#define COSMOLOGY_TOOLS_H_

#include "CosmologyToolsAPIMangling.h" // auto-generated for Fortran interface

#include "CosmologyToolsDefinitions.h"
#include "CosmologyToolsManager.h"
#include <mpi.h>

cosmotk::CosmologyToolsManager *CosmoToolsManager;

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initialize the cosmotools in-situ library.
 * @param comm user-supplied communicator where the analysis will take place.
 * @note This method is typically called once, in the beginning of the program.
 */
void cosmotools_initialize(MPI_Comm *comm);

/**
 * @brief Initialize the cosmotools in-situ library from a fortran program.
 * @param fcomm the user-supplied fortran communicator
 * @see cosmotools_initialize
 */
void cosmotools_fortran_initialize(MPI_Fint *fcomm);

/**
 * @brief Set the configuration file for the cosmotools library
 * @param configfile filename and path of the configuration file
 * @note This method is typically called once, in the beginning of the program.
 */
void cosmotools_set_analysis_config(char *configfile);


/**
 * @brief Set the domain parameters for the N-body
 * @param boxlength the length of the box domain
 * @param ghostoverlap the ghost zone overlap
 * @param NDIM the number of points along each dimension
 * @note This method is typically called once, in the beginning of the program.
 */
void cosmotools_set_domain_parameters(
      REAL boxlength, INTEGER ghostoverlap, INTEGER NDIM);

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
void cosmotools_set_particles(
        INTEGER *tstep, REAL *redshift,
        POSVEL_T *px, POSVEL_T *py, POSVEL_T *pz,
        POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz,
        REAL *mass, POTENTIAL_T *potential,
        ID_T *tags, MASK_T *mask, STATUS_T *status,
        ID_T *N);

/**
 * @brief Called at each time-step to perform analysis.
 * @note The analysis task(s) that will be executed at each time-step are
 * defined on the user-supplied configuration file. Each time coprocess is
 * called, rank 0 (w.r.t. to the given communicator) parses the configuration
 * file and constructs a CosmoToolsParameters object. Then CosmoToolsParameters
 * is distributed to all ranks.
 */
void cosmotools_coprocess();

/**
 * @brief Finalizes the cosmotools environment.
 */
void cosmotools_finalize();



#ifdef __cplusplus
}
#endif

#endif /* COSMOLOGY_TOOLS_H_ */
