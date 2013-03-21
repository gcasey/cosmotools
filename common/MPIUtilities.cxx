#include "MPIUtilities.h"

#include <cstdio>
#include <cstdarg>

namespace cosmotk
{

MPIUtilities::MPIUtilities()
{
  // TODO Auto-generated constructor stub

}

//------------------------------------------------------------------------------
MPIUtilities::~MPIUtilities()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void MPIUtilities::Printf(MPI_Comm comm, const char *fmt,...)
{
  // Sanity check
  assert("pre: MPI communicator should not be NULL!" && (comm != MPI_COMM_NULL) );
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
void MPIUtilities::SynchronizedPrintf(MPI_Comm comm, const char *fmt,...)
{
  // Sanity check
  assert("pre: MPI communicator should not be NULL!" && (comm != MPI_COMM_NULL) );
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
    printf("[%d]: ",rank);
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
    printf("[%d]: ",rank);
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
    printf("[%d]: ",rank);
    vprintf(fmt,argptr);
    fflush(stdout);
    va_end(argptr);

    // STEP 2: tell next process to print
    MPI_Isend(NULL,0,MPI_INT,rank+1,0,comm,&nullRequest);
    }
}

} /* namespace cosmotk */
