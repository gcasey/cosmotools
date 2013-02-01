/**
 * @brief Simple program developed to facilitate the design and development
 * for doing I/O on the MergerTree using DIY.
 */

// C++ includes
#include <iostream>

// MPI include
#include <mpi.h>

// DIY include
#include "diy.h"

// CosmologyTools includes
#include "DistributedHaloEvolutionTree.h"
#include "Halo.h"

//==============================================================================
// Global variables
//==============================================================================
int rank;
int size;
MPI_Comm comm;
cosmotk::DistributedHaloEvolutionTree HaloMergerTree;

DIYTree MTree;

std::string FileName = "MergerTree.dat";

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
void CreateMergerTree();
void WriteMergerTree();

void* create_datatype(void *block, int did, int lid, DIY_Datatype *dtype)
{
  cosmotk::DistributedHaloEvolutionTree::CreateDIYTreeType(&MTree,dtype);
  return DIY_BOTTOM;
}

//==============================================================================
// Program main
//==============================================================================
int main(int argc, char **argv)
{
  // STEP 0: Initialize MPI
  MPI_Init(&argc,&argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);
  PRINTLN("- Initialized MPI...[DONE]");

  // STEP 1: Initialize DIY
  DIY_Init(3,NULL,1,comm);
  PRINTLN("- Initialized DIY...[DONE]");

  // STEP 2: Decompose domain
  int tot_blocks = size; // total number of blocks is equal to number of ranks
  int nblocks    = 1;     // local blocks per process is 1
  int given[3]   = {0, 0, 0}; // no constraints on decomposition in {x, y, z}
  int ghost[6]   = {0, 0, 0, 0, 0, 0}; // -x, +x, -y, +y, -z, +z ghost
  DIY_Decompose(ROUND_ROBIN_ORDER,tot_blocks,&nblocks,1,ghost,given);
  PRINTLN("- Prescribed decomposition to DIY...[DONE]");

  // STEP 3: Create MergerTree at each process
  CreateMergerTree();
  PRINTLN("- Build fake merger-tree...[DONE]");

  // STEP 4: Get DIY Tree
  HaloMergerTree.GetDIYTree(&MTree);

  // STEP 5: Write Merger-Tree
  WriteMergerTree();
  PRINTLN("- Write merger-tree...[DONE]");

  // STEP 6: Finalize DIY
  DIY_Finalize();
  PRINTLN("- Finalized DIY...[DONE]");

  // STEP 7: Finalize MPI
  MPI_Finalize();
  return 0;
}

//------------------------------------------------------------------------------
void WriteMergerTree()
{
  int nblocks = 1;
  int **hdrs = new int*[nblocks];
  void **pmblocks;
  pmblocks = new void*[nblocks];
  for(int blkIdx=0; blkIdx < nblocks; ++blkIdx )
    {
    std::cerr << "Block: "    << blkIdx              << std::endl;
    std::cerr << "NumNodes: " << MTree.NumberOfNodes << std::endl;
    std::cerr << "NumEdges: " << MTree.NumberOfEdges << std::endl;

    hdrs[blkIdx]     = new int[DIY_MAX_HDR_ELEMENTS];
    hdrs[blkIdx][0]  = MTree.NumberOfNodes;
    hdrs[blkIdx][1]  = MTree.NumberOfEdges;
    pmblocks[blkIdx] = &MTree;
    }
  DIY_Write_open_all(0,const_cast<char*>(FileName.c_str()),0);
  DIY_Write_blocks_all(0,pmblocks, nblocks, hdrs, &create_datatype);
  DIY_Write_close_all(0);

  delete [] pmblocks;
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
