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
#include "GenericIO.h"

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

//==============================================================================
// Function prototypes
//==============================================================================
void CartCommInit(MPI_Comm comm);

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

  // STEP 1: Make a cartesian communicator
  CartCommInit(comm);
  PRINTLN("- Initialize cartesian communicator...[DONE]");

  // STEP 2: Read and print vars
  cosmotk::GenericIO reader(comm,file);
  reader.openAndReadHeader(false);

  std::vector<cosmotk::GenericIO::VariableInfo> variables;
  reader.getVariableInfo(variables);
  if( rank == 0 )
    {
    std::cout << "Number of variables: " << variables.size() << std::endl;
    for(int i=0; i < variables.size(); ++i)
      {
      std::cout << variables[i].Name << " " << variables[i].Size << std::endl;
      } // END for
    } // END if

  reader.close();


  // STEP 3: Finalize
  MPI_Finalize();
  return 0;
}

//------------------------------------------------------------------------------
void CartCommInit(MPI_Comm mycomm)
{
  int periodic[] = { 1, 1, 1 };
  int reorder = 0;
  int dims[3] = { 0, 0, 0 };

  // Get cartesian dimensions
  MPI_Dims_create(size,3,dims);

  // Create cartesian communicator
  MPI_Comm cartComm;
  MPI_Cart_create(mycomm,3,dims,periodic,reorder,&cartComm);

  // Reset global variables for rank and comm
  MPI_Comm_rank(cartComm,&rank);
  comm=cartComm;
}

