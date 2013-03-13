/**
 * @brief A simple program that opens a file written in the GenericIO format
 * and prints out the number of variables, their corresponding sizes and
 * names.
 */

// C++ includes
#include <iostream>
#include <vector>

// MPI include
#include <mpi.h>

// CosmologyTools includes
#include "GenericIOMPIReader.h"

//==============================================================================
// Global variables
//==============================================================================
int rank;
int size;
MPI_Comm comm = MPI_COMM_WORLD;

//==============================================================================
// Macros
//==============================================================================
#define PRINTLN( str ) {          \
  if(rank==0) {                    \
    std::cout << str << std::endl; \
    std::cout.flush();             \
  }                                \
  MPI_Barrier(comm);               \
}

#define PRINT( str ) {    \
  if(rank==0) {            \
    std::cout << str;      \
    std::cout.flush();     \
  }                        \
  MPI_Barrier(comm);       \
}

/**
 * @brief Program main
 * @param argc argument counter
 * @param argv argument vector
 * @return rc return code
 */
int main(int argc, char **argv)
{
  // STEP 0: Initialize
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  PRINTLN("- Initialize MPI...[DONE]");

  // STEP 1: Get file name to read
  std::string file=std::string(argv[1]);

  // STEP 2: Read and print vars
  cosmotk::GenericIOMPIReader reader;
  reader.SetCommunicator(comm);
  reader.SetFileName(file);

  reader.OpenAndReadHeader();
  std::cout << "NumElements: " << reader.GetNumberOfElements() << std::endl;
  std::cout.flush();

  if( rank == 0 )
    {
    for(int i=0; i < reader.GetNumberOfVariables(); ++i)
      {
      std::cout << reader.GetVariableName( i ) << "\t";
      std::cout.flush();
      }
    std::cout << std::endl;
    std::cout.flush();
    } // END if
  MPI_Barrier(comm);

  reader.Close();

  // STEP 3: Finalize
  MPI_Finalize();
  return 0;
}


