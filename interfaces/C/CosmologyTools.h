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
 * @brief Sets the memory layout that the interface will be using.
 * @param memoryLayout the memory layout
 */
void CosmologySetMemoryLayout(int* memoryLayout);

/**
 * @brief Sets the frequency at which the tracker will be invoked.
 * @param frequency the frequency for the halo tracker
 */
void CosmologySetTrackerFrequency(int *frequency);

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
