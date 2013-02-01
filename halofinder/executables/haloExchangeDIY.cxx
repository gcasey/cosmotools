/**
 * @brief Simple program developed to facilitate the design and development
 * for exchanging halos.
 */

// C++ includes
#include <iostream>  // For C++ I/O stream object
#include <fstream>   // For C++ file stream object
#include <sstream>   // For C++ string stream object
#include <vector>    // For STL vector

// MPI include
#include <mpi.h>

// DIY include
#include "diy.h"

// Cosmotools include
#include "CosmologyToolsMacros.h"
#include "Halo.h"

//==============================================================================
// Global variables
//==============================================================================
static const float BOXLENGTH=100.0;
static const int NDIM=100;

int rank;
int size;
MPI_Comm comm;

std::map<std::string,cosmotk::Halo> RcvHalos;

unsigned char neigh_dirs[] = {
  DIY_X0,                   DIY_X1,
  DIY_Y0,                   DIY_Y1,
  DIY_Z0,                   DIY_Z1,
  DIY_X0 | DIY_Y0,          DIY_X1 | DIY_Y1,
  DIY_X0 | DIY_Y1,          DIY_X1 | DIY_Y0,
  DIY_Y0 | DIY_Z0,          DIY_Y1 | DIY_Z1,
  DIY_Y0 | DIY_Z1,          DIY_Y1 | DIY_Z0,
  DIY_Z0 | DIY_X0,          DIY_Z1 | DIY_X1,
  DIY_Z0 | DIY_X1,          DIY_Z1 | DIY_X0,
  DIY_X0 | DIY_Y0 | DIY_Z0, DIY_X1 | DIY_Y1 | DIY_Z1,
  DIY_X0 | DIY_Y0 | DIY_Z1, DIY_X1 | DIY_Y1 | DIY_Z0,
  DIY_X0 | DIY_Y1 | DIY_Z0, DIY_X1 | DIY_Y0 | DIY_Z1,
  DIY_X0 | DIY_Y1 | DIY_Z1, DIY_X1 | DIY_Y0 | DIY_Z0,
};


//
// Neighbors are enumerated so that particles can be attached to the correct
// neighbor, but these pairs must be preserved for the ParticleExchange.
// Every processor should be able to send and receive on every iteration of
// the exchange, so if everyone sends RIGHT and receives LEFT it works
//
// Do not change this pairing order.
//
enum NEIGHBOR
{
  X0,                   // Left face
  X1,                   // Right face

  Y0,                   // Bottom face
  Y1,                   // Top face

  Z0,                   // Front face
  Z1,                   // Back face

  X0_Y0,                // Left   bottom edge
  X1_Y1,                // Right  top    edge

  X0_Y1,                // Left   top    edge
  X1_Y0,                // Right  bottom edge

  Y0_Z0,                // Bottom front  edge
  Y1_Z1,                // Top    back   edge

  Y0_Z1,                // Bottom back   edge
  Y1_Z0,                // Top    front  edge

  Z0_X0,                // Front  left   edge
  Z1_X1,                // Back   right  edge

  Z0_X1,                // Front  right  edge
  Z1_X0,                // Back   left   edge

  X0_Y0_Z0,             // Left  bottom front corner
  X1_Y1_Z1,             // Right top    back  corner

  X0_Y0_Z1,             // Left  bottom back  corner
  X1_Y1_Z0,             // Right top    front corner

  X0_Y1_Z0,             // Left  top    front corner
  X1_Y0_Z1,             // Right bottom back  corner

  X0_Y1_Z1,             // Left  top    back  corner
  X1_Y0_Z0              // Right bottom front corner
};

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
bool IsCartesianComm(MPI_Comm comm);

void MPICart2DIYDecomposition(MPI_Comm comm);
void GetBlockBounds(
    int decompSize[3], int pos[3], float min[3], float dx[3]);
void ComputeRankNeighbors(int pos[3], int neighbor[26]);
int GetRankByPosition(int i, int j, int k);

void GetProcessHalos(std::vector<cosmotk::Halo> &halos);

void ExchangeHalos(std::vector<cosmotk::Halo> &halos);
void ExchangeHaloInfo(std::vector<cosmotk::Halo> &halos);
void ExchangeHaloParticles(std::vector<cosmotk::Halo> &halos);

void WriteRcvHalos();


//==============================================================================
// Program main
//==============================================================================
int main(int argc, char **argv)
{
  // STEP 0: Initialize MPI
  MPI_Init(&argc, &argv);
  comm = MPI_COMM_WORLD;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &size);
  PRINTLN("- Initialized MPI...[DONE]");

  // STEP 1: Create cartesian communicator
  CartCommInit(comm);
  assert("pre: Could not create cartesian communicator" &&
          IsCartesianComm(comm));
  PRINTLN("- Defined cartesian communicator...[DONE]");


  // STEP 2: Initialize DIY
  DIY_Init(3, NULL, 1, comm);
  PRINTLN("- Initialized DIY...[DONE]");

  // STEP 3: Prescribe decomposition to DIY
  MPICart2DIYDecomposition(comm);
  PRINTLN("- Prescribed decomposition to DIY...[DONE]");

  // STEP 4: Get mock-halos in this process
  std::vector< cosmotk::Halo > halos;
  GetProcessHalos( halos );

  // STEP 5: Exchange halos
  PRINTLN("- Exchange halos...");
  ExchangeHalos( halos );
  PRINTLN("- Exchange halos...[DONE]");

  // STEP 6: Write receive halos
  PRINTLN("- Write Receive halos...");
  WriteRcvHalos();
  PRINTLN("- Write Receive halos...[DONE]");

  // STEP 6: Finalize MPI/DIY
  DIY_Finalize();
  PRINTLN("- Finalized DIY...[DONE]");
  MPI_Finalize();
  return 0;
}

//------------------------------------------------------------------------------
void WriteRcvHalos()
{
  std::ostringstream oss;
  oss << rank << "_rcv_halos.dat";

  std::ofstream ofs;
  ofs.open(oss.str().c_str());
  std::map<std::string,cosmotk::Halo>::iterator haloIter = RcvHalos.begin();
  for(; haloIter != RcvHalos.end(); ++haloIter )
    {
    ofs << "==============================\n";
    haloIter->second.Print(ofs);
    ofs << "\n";
    } // END for all halos
  ofs.close();
}

//------------------------------------------------------------------------------
void ExchangeHaloInfo(std::vector<cosmotk::Halo> &halos)
{
  DIYHaloItem haloInfo;
  for( int hidx=0; hidx < halos.size(); ++hidx )
    {
    halos[ hidx ].GetDIYHaloItem( &haloInfo );
    DIY_Enqueue_item_all(
        0, 0, (void *)&haloInfo, NULL, sizeof(DIYHaloItem), NULL);
    } // END for all halos
  PRINTLN("Enqueued halos!");

  int nblocks = 1;
  void ***rcvHalos      = new void**[nblocks];
  int *numHalosReceived = new int[nblocks];
  DIY_Exchange_neighbors(
      0, rcvHalos,numHalosReceived,1.0,&cosmotk::Halo::CreateDIYHaloType);

  DIYHaloItem *rcvHaloItem = NULL;
  for(int i=0; i < nblocks; ++i)
    {
    for( int j=0; j < numHalosReceived[i]; ++j)
      {
      rcvHaloItem = (struct DIYHaloItem *)rcvHalos[i][j];
      cosmotk::Halo h(rcvHaloItem);
      RcvHalos[ h.GetHashCode() ] = h;
      } // END for all received halos of this block
    } // END for all blocks

  // clean up
  DIY_Flush_neighbors(
      0, rcvHalos,numHalosReceived,&cosmotk::Halo::CreateDIYHaloType);
  delete [] numHalosReceived;
}

//------------------------------------------------------------------------------
void ExchangeHaloParticles(std::vector<cosmotk::Halo> &halos)
{
  std::vector<DIYHaloParticleItem> haloParticles;
  for(int hidx=0; hidx < halos.size(); ++hidx)
    {
    haloParticles.resize(0);
    halos[hidx].GetDIYHaloParticleItemsVector(haloParticles);
    for( int pIdx=0; pIdx < haloParticles.size(); ++pIdx )
      {
      DIY_Enqueue_item_all(
          0,0,
          (void *)&haloParticles[pIdx],
          NULL,
          sizeof(DIYHaloParticleItem),
          NULL);
      } // END for all particles within the halo
    } // END for all halos

  int nblocks = 1;
  void ***rcvHalos = new void**[nblocks];
  int *numHalosReceived = new int[nblocks];
  DIY_Exchange_neighbors(
      0,rcvHalos,numHalosReceived,1.0,&cosmotk::Halo::CreateDIYHaloParticleType);

  DIYHaloParticleItem *haloParticle = NULL;
  for(int i=0; i < nblocks; ++i)
    {
    for(int j=0; j < numHalosReceived[i]; ++j)
      {
      haloParticle = (struct DIYHaloParticleItem*)rcvHalos[i][j];
      std::string hashCode =
          cosmotk::Halo::GetHashCodeForHalo(
              haloParticle->Tag,haloParticle->TimeStep);
      assert(RcvHalos.find(hashCode)!=RcvHalos.end());
      RcvHalos[hashCode].ParticleIds.insert(haloParticle->HaloParticleID);
      } // END for all received halos of this block
    } // END for all blocks

  // clean up
  DIY_Flush_neighbors(
      0,rcvHalos,numHalosReceived,&cosmotk::Halo::CreateDIYHaloParticleType);
  delete [] numHalosReceived;
}

//------------------------------------------------------------------------------
void ExchangeHalos(std::vector< cosmotk::Halo > &halos )
{
  PRINT("- Exchange halo information...");
  ExchangeHaloInfo( halos );
  PRINTLN("[DONE]");

  PRINT("- Exchange halo particles...");
  ExchangeHaloParticles( halos );
  PRINTLN("[DONE]");
}

//------------------------------------------------------------------------------
void GetProcessHalos(std::vector< cosmotk::Halo > &halos)
{
  int numHalosToCreate = 2;
  int numParticles     = 5;

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
}

//------------------------------------------------------------------------------
int GetRankByPosition(int i, int j, int k)
{
  int ijk[3];
  ijk[0]=i; ijk[1]=j; ijk[2]=k;
  int rnk;
  MPI_Cart_rank(comm,ijk,&rnk);
  return( rnk );
}

//------------------------------------------------------------------------------
void GetBlockBounds(
    int decompSize[3], int pos[3], float min[3], float dx[3])
{
  float m_grid2phys_pos = BOXLENGTH/static_cast<float>(NDIM);

  for(int i=0; i < 3; ++i)
    {
    dx[i] = BOXLENGTH/static_cast<float>(decompSize[i]);
    float m_corner_grid_alive = pos[i]*(NDIM/decompSize[i]);
    min[i] = m_corner_grid_alive*m_grid2phys_pos;
    } // END for
}

//------------------------------------------------------------------------------
void ComputeRankNeighbors(int pos[3], int neighbor[26])
{
 assert("pre: communicator is not cartesian!" && IsCartesianComm(comm));

 int xpos = pos[0];
 int ypos = pos[1];
 int zpos = pos[2];

 // Face neighbors
 neighbor[X0] = GetRankByPosition(xpos-1, ypos, zpos);
 neighbor[X1] = GetRankByPosition(xpos+1, ypos, zpos);
 neighbor[Y0] = GetRankByPosition(xpos, ypos-1, zpos);
 neighbor[Y1] = GetRankByPosition(xpos, ypos+1, zpos);
 neighbor[Z0] = GetRankByPosition(xpos, ypos, zpos-1);
 neighbor[Z1] = GetRankByPosition(xpos, ypos, zpos+1);

 // Edge neighbors
 neighbor[X0_Y0] = GetRankByPosition(xpos-1, ypos-1, zpos);
 neighbor[X0_Y1] = GetRankByPosition(xpos-1, ypos+1, zpos);
 neighbor[X1_Y0] = GetRankByPosition(xpos+1, ypos-1, zpos);
 neighbor[X1_Y1] = GetRankByPosition(xpos+1, ypos+1, zpos);

 neighbor[Y0_Z0] = GetRankByPosition(xpos, ypos-1, zpos-1);
 neighbor[Y0_Z1] = GetRankByPosition(xpos, ypos-1, zpos+1);
 neighbor[Y1_Z0] = GetRankByPosition(xpos, ypos+1, zpos-1);
 neighbor[Y1_Z1] = GetRankByPosition(xpos, ypos+1, zpos+1);

 neighbor[Z0_X0] = GetRankByPosition(xpos-1, ypos, zpos-1);
 neighbor[Z0_X1] = GetRankByPosition(xpos+1, ypos, zpos-1);
 neighbor[Z1_X0] = GetRankByPosition(xpos-1, ypos, zpos+1);
 neighbor[Z1_X1] = GetRankByPosition(xpos+1, ypos, zpos+1);

 // Corner neighbors
 neighbor[X0_Y0_Z0] = GetRankByPosition(xpos-1, ypos-1, zpos-1);
 neighbor[X1_Y0_Z0] = GetRankByPosition(xpos+1, ypos-1, zpos-1);
 neighbor[X0_Y1_Z0] = GetRankByPosition(xpos-1, ypos+1, zpos-1);
 neighbor[X1_Y1_Z0] = GetRankByPosition(xpos+1, ypos+1, zpos-1);
 neighbor[X0_Y0_Z1] = GetRankByPosition(xpos-1, ypos-1, zpos+1);
 neighbor[X1_Y0_Z1] = GetRankByPosition(xpos+1, ypos-1, zpos+1);
 neighbor[X0_Y1_Z1] = GetRankByPosition(xpos-1, ypos+1, zpos+1);
 neighbor[X1_Y1_Z1] = GetRankByPosition(xpos+1, ypos+1, zpos+1);
}

//------------------------------------------------------------------------------
void MPICart2DIYDecomposition(MPI_Comm mycomm)
{
  assert("pre: MPI communicator is NULL!" &&
          (mycomm != MPI_COMM_NULL) );
  assert("pre: communicator is not cartesian!" &&
          IsCartesianComm(mycomm));

  // STEP 0: We are dealing with a 3-D domain and with structured topology, so,
  // each block has 26 neighbors.
  const int DIMENSION        = 3;
  const int NUM_OF_NEIGHBORS = 26;

  // STEP 1: Get cartesian topology
  int myPosition[3];  // the position of this rank
  int decomp_size[3]; // the dimensions of the cartesian topology
  int periodicity[3]; // periodicity for each direction
  MPI_Cart_get(
     mycomm,DIMENSION,decomp_size,periodicity,myPosition);

  // STEP 2: Compute number of blocks per process and total number of blocks
  // Currently, numBlocksPerProcess = 1 and totalNumberOfBlocks is given by
  // numBlocksPerProcess * numRanks.
  int numBlocksPerProcess = 1;
  int totalBlocks = 0;
  int numRanks = size;
  totalBlocks = numBlocksPerProcess * numRanks;

  // STEP 3: Get the Block bounds
  bb_t bb;
  float min[DIMENSION], dx[DIMENSION];
  GetBlockBounds(decomp_size,myPosition,min,dx);
  for( int i=0; i < DIMENSION; ++i )
   {
   bb.min[i] = min[i];
   bb.max[i] = min[i] + dx[i];
   }

  // STEP 4: data overall extents
  // assume all blocks are same size (as mine)
  // assume 0,0,0 is the overall data minimum corner
  float data_mins[DIMENSION], data_maxs[DIMENSION];
  for (int i = 0; i < DIMENSION; i++) {
    data_mins[i] = 0.0;
    data_maxs[i] = data_mins[i] + dx[i] * decomp_size[i];
  }

  // STEP 5: Get neighbors
  int neigh_gids[NUM_OF_NEIGHBORS];
  ComputeRankNeighbors(myPosition,neigh_gids);
  gb_t **neighbors = new gb_t*[1];
  neighbors[0] = new gb_t[NUM_OF_NEIGHBORS];
  int num_neighbors[1] = { NUM_OF_NEIGHBORS };
  for (int i = 0; i < NUM_OF_NEIGHBORS; i++)
   {
   neighbors[0][i].gid       = neigh_gids[i];
   neighbors[0][i].proc      = neigh_gids[i];
   neighbors[0][i].neigh_dir = neigh_dirs[i];
   } // END for all neighbors

  // STEP 6: gids are trivial, only one block, my MPI rank
  int gids[1];
  gids[0] = rank;

  // STEP 7: Get wrap, i.e., domain is all-periodic
  int wrap = 1;

  // STEP 8: Describe the already decomposed domain in DIY
  DIY_Decomposed(numBlocksPerProcess,gids,&bb,
                 NULL,NULL,NULL,NULL,
                 neighbors,num_neighbors,wrap);
  delete [] neighbors[0];
  delete [] neighbors;
}

//------------------------------------------------------------------------------
bool IsCartesianComm(MPI_Comm mycomm)
{
  assert("pre: MPI communicator is NULL!" &&
          (mycomm != MPI_COMM_NULL) );

   int topology = 0;
   MPI_Topo_test(mycomm,&topology);
   return( (topology==MPI_CART) );
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
