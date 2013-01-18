/**
 * @brief Simple program developed to facilitate the design and development
 * for doing I/O on the MergerTree using the GenericIO infrastructure.
 */

// C++ includes
#include <fstream>
#include <iostream>
#include <sstream>

// MPI include
#include <mpi.h>

#include "diy.h" // Needed b/c the cosmotools framework has native DIY support

// CosmologyTools includes
#include "DistributedHaloEvolutionTree.h"
#include "GenericIO.h"
#include "Halo.h"

//==============================================================================
// Global variables
//==============================================================================
int rank;
int size;
MPI_Comm comm = MPI_COMM_WORLD;
std::string FileName = "MergerTree.dat";
cosmotk::DistributedHaloEvolutionTree HaloMergerTree;

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
void ReadMergerTree();


//==============================================================================
// Program main
//==============================================================================
int main(int argc, char **argv)
{
  // STEP 0: Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  PRINTLN("- Initialize MPI...[DONE]");

  // STEP 1: Initialize cartesian communicator
  CartCommInit(comm);
  PRINTLN("- Initialize cartesian communicator...[DONE]");

  // STEP 4: Read Merger-Tree
  ReadMergerTree();
  PRINTLN("- Read merger-tree...[DONE]");

  std::ostringstream oss;
  oss << "tree-" << rank << ".txt";

  std::ofstream ofs;
  ofs.open(oss.str().c_str());
  ofs << "=======================" << std::endl;
  ofs << HaloMergerTree.ToString();
  ofs << std::endl;
  ofs << "=======================" << std::endl;
  ofs.close();

  // STEP 5: Finalize MPI
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

//------------------------------------------------------------------------------
void ReadMergerTree()
{
  HaloMergerTree.SetFileName("mtree");
  HaloMergerTree.SetCommunicator(comm);
  HaloMergerTree.Clear();
  HaloMergerTree.ReadTree();
}





