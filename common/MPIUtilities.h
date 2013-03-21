/**
 * @brief A singleton class that implements generic functions for MPI apps/libs.
 */
#ifndef MPIUTILITIES_H_
#define MPIUTILITIES_H_

#include "CosmologyToolsMacros.h"

#include <mpi.h>

namespace cosmotk
{

class MPIUtilities
{
public:

  /**
   * @brief Rank 0 prints the formatted message to stdout
   * @param comm the communicator used by the parallel application
   * @param fmt the formatted message to print
   * @note This is a collective call, all processes must call this method.
   * @pre comm != MPI_COMM_NULL
   * @pre fmt != NULL
   */
  static void Printf(MPI_Comm comm, const char *fmt,...);

  /**
   * @brief Each process prints the formatted message in stdout in rank order.
   * That is, process "i" prints its message right after process "i-1".
   * @param comm the communicator used by the parallel application
   * @param fmt the formatted message to print
   * @note This is a collective call, all processes must call this method.
   * @pre comm != MPI_COMM_NULL
   * @pre fmt != NULL
   */
  static void SynchronizedPrintf(MPI_Comm comm, const char *fmt,...);

protected:
  MPIUtilities();
  virtual ~MPIUtilities();

private:
  DISABLE_COPY_AND_ASSIGNMENT(MPIUtilities);
};

} /* namespace cosmotk */
#endif /* MPIUTILITIES_H_ */
