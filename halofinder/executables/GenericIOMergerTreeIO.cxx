/**
 * @brief Simple program developed to facilitate the design and development
 * for doing I/O on the MergerTree using the GenericIO infrastructure.
 */

// C++ includes
#include <iostream>
#include <fstream>
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
void CreateMergerTree();
void WriteMergerTree();
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

  DIY_Init(3,NULL,1,comm);

  int tot_blocks = size; // total number of blocks is equal to number of ranks
  int nblocks    = 1;     // local blocks per process is 1
  int given[3]   = {0, 0, 0}; // no constraints on decomposition in {x, y, z}
  int ghost[6]   = {0, 0, 0, 0, 0, 0}; // -x, +x, -y, +y, -z, +z ghost
  DIY_Decompose(ROUND_ROBIN_ORDER,tot_blocks,&nblocks,1,ghost,given);

  // STEP 2: Create MergerTree at each process
  CreateMergerTree();
  PRINTLN("- Build fake merger-tree...[DONE]");

  std::ostringstream oss;
  oss << "expected-tree-" << rank << ".txt";

  std::ofstream ofs;
  ofs.open(oss.str().c_str());
  ofs << "=======================" << std::endl;
  ofs << HaloMergerTree.ToString();
  ofs << std::endl;
  ofs << "=======================" << std::endl;
  ofs.close();

  // STEP 3: Write Merger-Tree
  WriteMergerTree();
  PRINTLN("- Write merger-tree...[DONE]");

  // STEP 4: Read Merger-Tree
  ReadMergerTree();
  PRINTLN("- Read merger-tree...[DONE]");

  DIY_Finalize();

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
void WriteMergerTree()
{
  HaloMergerTree.SetFileName("mtree");
  HaloMergerTree.WriteTree();
}

//------------------------------------------------------------------------------
void ReadMergerTree()
{
  HaloMergerTree.Clear();
  HaloMergerTree.SetFileName("mtree");
  HaloMergerTree.ReadTree();
}

//------------------------------------------------------------------------------
void CreateMergerTree()
{
  HaloMergerTree.SetCommunicator(comm);

  int numHalosToCreate = 2;
  int numParticles     = 5;

  std::vector< cosmotk::Halo > halos;
  halos.resize(numHalosToCreate);
  for( int haloIdx=0; haloIdx < numHalosToCreate; ++haloIdx )
    {
    halos[haloIdx].TimeStep = 200;
    halos[haloIdx].Redshift = 20;
    halos[haloIdx].Tag      = rank*10+haloIdx;

    for( int pIdx=0; pIdx < numParticles; ++pIdx )
      {
      ID_T particleID = halos[haloIdx].Tag + pIdx;
      halos[haloIdx].ParticleIds.insert(particleID);
      } // END for all halo particles
    } // END for all halos

  HaloMergerTree.AppendNodes( &halos[0], halos.size() );

  for(int haloIdx=0; haloIdx < halos.size()-1; ++haloIdx)
    {
    HaloMergerTree.CreateEdge(
        halos[haloIdx].GetHashCode(),
        halos[haloIdx+1].GetHashCode() );
    }

  HaloMergerTree.RelabelTreeNodes();
}

