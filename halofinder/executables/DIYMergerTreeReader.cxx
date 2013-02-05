/**
 * @brief Simple program developed to facilitate the design and development
 * for doing I/O on the MergerTree using DIY.
 */

// C++ includes
#include <iostream>
#include <fstream>

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
void ReadData();

void* create_datatype(int did, int lid, int *hdr, DIY_Datatype *dtype)
{
  std::cerr << "Header: " << hdr[0] << " " << hdr[1] << std::endl;
  MTree.NumberOfNodes = hdr[0];
  MTree.NumberOfEdges = hdr[1];
  MTree.TreeNodes = new DIYTreeNodeType[MTree.NumberOfNodes];
  MTree.TreeEdges = new DIYTreeEdgeType[MTree.NumberOfEdges];
  cosmotk::DistributedHaloEvolutionTree::CreateDIYTreeType(&MTree,dtype);
  return &MTree;
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

  // STEP 3: Read the data
  ReadData();

  // STEP 4: Finalize DIY
  DIY_Finalize();
  PRINTLN("- Finalized DIY...[DONE]");

  // STEP 5: Finalize MPI
  MPI_Finalize();
  return 0;
}

//==============================================================================
void ReadData()
{
  int swap     = 0;
  int compress = 0;
  int numblocks =
   DIY_Read_open_all(0,const_cast<char *>(FileName.c_str()),swap,compress);
  PRINTLN("- Open output file...[DONE]");

  void **pmblocks;
  int **hdrs = new int*[numblocks];
  for(int i=0; i < numblocks; ++i )
    {
    hdrs[i] = new int[DIY_MAX_HDR_ELEMENTS];
    }

  PRINT("- Read blocks...");
  DIY_Read_blocks_all(0,&pmblocks,hdrs,&create_datatype);
  PRINTLN("[DONE]");

  PRINT("- Close file...");
  DIY_Read_close_all(0);
  PRINTLN("[DONE]");

  std::ofstream ofs;
  ofs.open("MergerTreeAscii.dat");
  ofs << "NUMBER OF BLOCKS " << numblocks           << std::endl;
  ofs << "NUMBER OF NODES "  << MTree.NumberOfNodes << std::endl;
  for(int nodeIdx=0; nodeIdx < MTree.NumberOfNodes; ++nodeIdx)
    {
    ofs << MTree.TreeNodes[nodeIdx].UniqueNodeID    << "\t";
    ofs << MTree.TreeNodes[nodeIdx].HaloTag         << "\t";
    ofs << MTree.TreeNodes[nodeIdx].TimeStep        << "\t";
    ofs << MTree.TreeNodes[nodeIdx].HaloCenter[0]   << ",";
    ofs << MTree.TreeNodes[nodeIdx].HaloCenter[1]   << ",";
    ofs << MTree.TreeNodes[nodeIdx].HaloCenter[2]   << "\t";
    ofs << MTree.TreeNodes[nodeIdx].HaloVelocity[0] << ",";
    ofs << MTree.TreeNodes[nodeIdx].HaloVelocity[1] << ",";
    ofs << MTree.TreeNodes[nodeIdx].HaloVelocity[2] << "\t";
    ofs << std::endl;
    } // END for all nodes
  ofs << std::endl << std::endl;
  ofs << "NUMBER OF EDGES "  << MTree.NumberOfEdges << std::endl;
  for(int edgeIdx=0; edgeIdx < MTree.NumberOfEdges; ++edgeIdx)
    {
    ofs << MTree.TreeEdges[edgeIdx].EndNodes[0] << ",";
    ofs << MTree.TreeEdges[edgeIdx].EndNodes[1] << "\t";
    ofs << MTree.TreeEdges[edgeIdx].EdgeWeight  << "\t";
    ofs << MTree.TreeEdges[edgeIdx].EdgeEvent   << "\t";
    ofs << std::endl;
    } // END for all edges
  ofs.close();
}
