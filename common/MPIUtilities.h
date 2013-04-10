/**
 * @brief A singleton class that implements generic functions for MPI apps/libs.
 */
#ifndef MPIUTILITIES_H_
#define MPIUTILITIES_H_

#include "CosmologyToolsMacros.h"

#include <mpi.h>

namespace cosmotk {
  namespace MPIUtilities {

/**
 * @brief Given a distributed set of items, wherein each process, P_i, has
 * a variable number of items, numItems, this method computes a range for
 * each process that corresponds to a global identifier for each item.
 * @param comm the underlying MPI communicator to use (in).
 * @param numItems the number of items owned by this process (in).
 * @param range the range assigned to this process (in/out).
 * @pre comm != MPI_COMM_NULL
 */
void GetProcessRange(
        MPI_Comm comm, ID_T numItems, ID_T range[2]);

/**
 * @brief Rank 0 prints the formatted message to stdout
 * @param comm the communicator used by the parallel application
 * @param fmt the formatted message to print
 * @note This is a collective call, all processes must call this method.
 * @pre comm != MPI_COMM_NULL
 * @pre fmt != NULL
 */
void Printf(MPI_Comm comm, const char *fmt,...);

/**
 * @brief Each process prints the formatted message in stdout in rank order.
 * That is, process "i" prints its message right after process "i-1".
 * @param comm the communicator used by the parallel application
 * @param fmt the formatted message to print
 * @note This is a collective call, all processes must call this method.
 * @pre comm != MPI_COMM_NULL
 * @pre fmt != NULL
 */
void SynchronizedPrintf(MPI_Comm comm, const char *fmt,...);


  } /* namespace MPIUtilities */
} /* namespace cosmotk */
#endif /* MPIUTILITIES_H_ */
