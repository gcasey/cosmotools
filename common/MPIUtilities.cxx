#include "MPIUtilities.h"

#include <cassert>
#include <cstdarg>
#include <cstdio>

namespace cosmotk {
  namespace MPIUtilities {


//------------------------------------------------------------------------------
void GetProcessRange(
      MPI_Comm comm, ID_T numItems, ID_T range[2])
{
  // Sanity check
  assert( "pre: MPI communicator should not be NULL!" &&
          (comm != MPI_COMM_NULL) );

  // STEP 0: Initialize
  range[0] = range[1] = 0;

  // STEP 1: Compute prefix sum
  ID_T prefixSum = -1;
  MPI_Scan(&numItems,&prefixSum,1,MPI_ID_T,MPI_SUM,comm);

  // STEP 2: Get the range in this process
  range[0] = (prefixSum-1)-numItems+1;
  range[1] = prefixSum-1;
}

//------------------------------------------------------------------------------
void Printf(MPI_Comm comm, const char *fmt,...)
{
  // Sanity check
  assert("pre: MPI communicator should not be NULL!" &&
          (comm != MPI_COMM_NULL) );
  assert("pre: NULL formatted message!" && (fmt != NULL) );

  int rank = -1;
  MPI_Comm_rank(comm,&rank);
  if( rank==0 )
    {
    va_list argptr;
    va_start(argptr,fmt);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);
    }
  MPI_Barrier(comm);
}

//------------------------------------------------------------------------------
void SynchronizedPrintf(MPI_Comm comm, const char *fmt,...)
{
  // Sanity check
  assert("pre: MPI communicator should not be NULL!" &&
          (comm != MPI_COMM_NULL) );
  assert("pre: NULL formatted message!" && (fmt != NULL) );

  int rank     = -1;
  int numRanks = -1;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&numRanks);
  MPI_Request nullRequest = MPI_REQUEST_NULL;

  if( rank == 0 )
    {
    // This is the first rank

    // STEP 0: Print message
    va_list argptr;
    va_start(argptr,fmt);
    printf("# [%d]: ",rank);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);

    // STEP 1: Tell next process to print, if any
    if( numRanks > 1)
      {
      MPI_Isend(NULL,0,MPI_INT,rank+1,0,comm,&nullRequest);
      }
    }
  else if( rank == numRanks-1 )
    {
    // Last rank

    // STEP 0: Block until previous process completes
    MPI_Recv(NULL,0,MPI_INT,rank-1,MPI_ANY_TAG,comm,MPI_STATUS_IGNORE);

    // STEP 1: Print message
    va_list argptr;
    va_start(argptr,fmt);
    printf("# [%d]: ",rank);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);
    }
  else
    {
    // STEP 0: Block until previous process completes
    MPI_Recv(NULL,0,MPI_INT,rank-1,MPI_ANY_TAG,comm,MPI_STATUS_IGNORE);

    // STEP 1: Print message
    va_list argptr;
    va_start(argptr,fmt);
    printf("# [%d]: ",rank);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);

    // STEP 2: tell next process to print
    MPI_Isend(NULL,0,MPI_INT,rank+1,0,comm,&nullRequest);
    }
}

  } /*namespace MPIUtilities */
} /* namespace cosmotk */
