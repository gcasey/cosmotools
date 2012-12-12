#include "TessVoidFinderAnalysisTool.h"

#include "diy.h"
#include "tess.h"
#include "voronoi.h"

#include <cassert>
#include <iostream>
#include <sstream>

#include <mpi.h>

namespace cosmotk
{

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


TessVoidFinderAnalysisTool::TessVoidFinderAnalysisTool()
{
  this->Name         = "TESS";
  this->Communicator = MPI_COMM_NULL;
  this->MinVol = .0001;
  this->MaxVol = -1.0;
  this->GhostFactor = 5.0;
  this->CellSize  = 1.0;
  this->Initialized = false;
  this->TimeStatistics = new double[MAX_TIMES];
}

//------------------------------------------------------------------------------
TessVoidFinderAnalysisTool::~TessVoidFinderAnalysisTool()
{
  tess_finalize();
  if( TimeStatistics != NULL )
    {
    delete [] this->TimeStatistics;
    }
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::ParseParameters()
{
  // STEP 0: parse basic parameters, as defined by super-class
  this->ParseBasicParameters();

  // STEP 1: parse tess parameters
  this->GhostFactor =
      static_cast<REAL>(this->GetDoubleParameter("GHOST_ZONE_SIZE"));
  this->MinVol =
      static_cast<REAL>(this->GetDoubleParameter("MIN_VOL_THRESHOLD"));
  this->MaxVol =
      static_cast<REAL>(this->GetDoubleParameter("MAX_VOL_THRESHOLD"));
  this->CellSize =
      static_cast<REAL>(this->GetDoubleParameter("CELL_SIZE"));
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::Execute(SimulationParticles *particles)
{
  assert("pre: input particles are NULL!" && (particles != NULL));
  assert("pre: MPI communicator is NULL!" &&
         (this->Communicator != MPI_COMM_NULL) );

  // STEP 0: short-circuit here
  if( particles->NumParticles == 0 )
    {
    return;
    }

  // STEP 1: parse the analysis tool parameters
  this->ParseParameters();

  // STEP 2: Package particle positions for tess
  int num_points[1];
  num_points[0] = particles->NumParticles;
  float **positions = NULL;
  this->PackageParticlePositions(particles, &positions);

  // STEP 3: initialize tess, i.e., call tess_init()
  if( !this->Initialized )
    {
    this->InitializeTess(particles);
    }

  // STEP 4: execute tess on the given particle dataset, i.e., call tess()
  std::ostringstream oss;
  oss.str(""); oss.clear();
  oss << this->OutputFile << "-" << particles->TimeStep << ".out";
  tess(positions,num_points,const_cast<char*>(oss.str().c_str()));

  // STEP 5: Clear packaged particles
  this->ClearParticlePositions(positions);
  assert("post: 2-D positions array must be NULL!" && (positions==NULL));
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::WriteOutput()
{
  // NOTE: tess combines I/O in TessVoidFinderAnalysisTool::Execute
  return;
}

//------------------------------------------------------------------------------
std::string TessVoidFinderAnalysisTool::GetInformation()
{
  return(this->GetBasicInformation());
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::InitializeTess(
        SimulationParticles *particles)
{
  assert("pre: input particles are NULL!" && (particles != NULL) );
  assert("pre: MPI communicator is NULL!" &&
         (this->Communicator != MPI_COMM_NULL) );

  // STEP 0: We are dealing with a 3-D domain and with structured topology, so,
  // each block has 26 neighbors.
  const int DIMENSION        = 3;
  const int NUM_OF_NEIGHBORS = 26;

  // STEP 1: Ensure we are dealing with a cartesian communicator
  int topology = 0;
  MPI_Topo_test(this->Communicator,&topology);
  assert("pre: communicator is not cartesian!" && topology==MPI_CART);
  if( topology != MPI_CART )
    {
    std::cerr << "ERROR: tess is expecting a cartesian communicator!\n";
    MPI_Abort(this->Communicator,-1);
    }

  // STEP 2: Get cartesian topology
  int myPosition[DIMENSION];  // the position of this rank
  int decomp_size[DIMENSION]; // the dimensions of the cartesian topology
  int periodicity[DIMENSION]; // periodicity for each direction
  MPI_Cart_get(
      this->Communicator,DIMENSION,decomp_size,periodicity,myPosition);

  // STEP 3: Compute number of blocks per process and total number of blocks
  // Currently, numBlocksPerProcess = 1 and totalNumberOfBlocks is given by
  // numBlocksPerProcess * numRanks.
  int numBlocksPerProcess = 1;
  int totalBlocks = 0;
  int numRanks = 0;
  int rank = 0;
  MPI_Comm_rank(this->Communicator,&rank);
  MPI_Comm_size(this->Communicator,&numRanks);
  totalBlocks = numBlocksPerProcess * numRanks;

  // STEP 4: Get the Block bounds
  bb_t bb;
  float min[DIMENSION], size[DIMENSION];
  this->GetBlockBounds(decomp_size,myPosition,min,size);
  for( int i=0; i < DIMENSION; ++i )
    {
    bb.min[i] = min[i];
    bb.max[i] = min[i] + size[i];
    }

  // STEP 5: data overall extents
  // assume all blocks are same size (as mine)
  // assume 0,0,0 is the overall data minimum corner
  float data_mins[DIMENSION], data_maxs[DIMENSION];
  for (int i = 0; i < DIMENSION; i++) {
    data_mins[i] = 0.0;
    data_maxs[i] = data_mins[i] + size[i] * decomp_size[i];
  }

  // STEP 6: Get neighbors
  int neigh_gids[NUM_OF_NEIGHBORS];
  this->ComputeRankNeighbors(myPosition,neigh_gids);
  gb_t **neighbors = new gb_t*[1];
  neighbors[0] = new gb_t[NUM_OF_NEIGHBORS];
  int num_neighbors[1] = { NUM_OF_NEIGHBORS };
  for (int i = 0; i < NUM_OF_NEIGHBORS; i++)
    {
    neighbors[0][i].gid       = neigh_gids[i];
    neighbors[0][i].proc      = neigh_gids[i];
    neighbors[0][i].neigh_dir = neigh_dirs[i];
    } // END for all neighbors

  // STEP 7: gids are trivial, only one block, my MPI rank
  int gids[1];
  gids[0] = rank;

  // STEP 8: Get wrap
  int wrap = (this->Periodic)? 1:0;

  // STEP 9: call tess_init
  tess_init(
      numBlocksPerProcess,gids,&bb,neighbors,num_neighbors,
      this->CellSize,this->GhostFactor,data_mins,data_maxs,wrap,
      this->MinVol,this->MaxVol,this->Communicator,
      this->TimeStatistics);

  // STEP 10: Return all dynamically allocated memory
  delete[] neighbors[0];
  delete[] neighbors;

  // STEP 11: Set initialized to true
  this->Initialized = true;
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::ComputeRankNeighbors(
    int pos[3], int neighbor[26])
{
  int xpos = pos[0];
  int ypos = pos[1];
  int zpos = pos[2];

  // Face neighbors
  neighbor[X0] = this->GetRankByPosition(xpos-1, ypos, zpos);
  neighbor[X1] = this->GetRankByPosition(xpos+1, ypos, zpos);
  neighbor[Y0] = this->GetRankByPosition(xpos, ypos-1, zpos);
  neighbor[Y1] = this->GetRankByPosition(xpos, ypos+1, zpos);
  neighbor[Z0] = this->GetRankByPosition(xpos, ypos, zpos-1);
  neighbor[Z1] = this->GetRankByPosition(xpos, ypos, zpos+1);

  // Edge neighbors
  neighbor[X0_Y0] = this->GetRankByPosition(xpos-1, ypos-1, zpos);
  neighbor[X0_Y1] = this->GetRankByPosition(xpos-1, ypos+1, zpos);
  neighbor[X1_Y0] = this->GetRankByPosition(xpos+1, ypos-1, zpos);
  neighbor[X1_Y1] = this->GetRankByPosition(xpos+1, ypos+1, zpos);

  neighbor[Y0_Z0] = this->GetRankByPosition(xpos, ypos-1, zpos-1);
  neighbor[Y0_Z1] = this->GetRankByPosition(xpos, ypos-1, zpos+1);
  neighbor[Y1_Z0] = this->GetRankByPosition(xpos, ypos+1, zpos-1);
  neighbor[Y1_Z1] = this->GetRankByPosition(xpos, ypos+1, zpos+1);

  neighbor[Z0_X0] = this->GetRankByPosition(xpos-1, ypos, zpos-1);
  neighbor[Z0_X1] = this->GetRankByPosition(xpos+1, ypos, zpos-1);
  neighbor[Z1_X0] = this->GetRankByPosition(xpos-1, ypos, zpos+1);
  neighbor[Z1_X1] = this->GetRankByPosition(xpos+1, ypos, zpos+1);

  // Corner neighbors
  neighbor[X0_Y0_Z0] = this->GetRankByPosition(xpos-1, ypos-1, zpos-1);
  neighbor[X1_Y0_Z0] = this->GetRankByPosition(xpos+1, ypos-1, zpos-1);
  neighbor[X0_Y1_Z0] = this->GetRankByPosition(xpos-1, ypos+1, zpos-1);
  neighbor[X1_Y1_Z0] = this->GetRankByPosition(xpos+1, ypos+1, zpos-1);
  neighbor[X0_Y0_Z1] = this->GetRankByPosition(xpos-1, ypos-1, zpos+1);
  neighbor[X1_Y0_Z1] = this->GetRankByPosition(xpos+1, ypos-1, zpos+1);
  neighbor[X0_Y1_Z1] = this->GetRankByPosition(xpos-1, ypos+1, zpos+1);
  neighbor[X1_Y1_Z1] = this->GetRankByPosition(xpos+1, ypos+1, zpos+1);
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::GetBlockBounds(
      int decompSize[3], int pos[3], float min[3], float size[3])
{
  // NOTE: This method essentially computes the following from HACC
  // Domain::rL_local_alive(size);
  // Domain::corner_phys_alive(min);

  float m_grid2phys_pos = this->BoxLength/static_cast<float>(this->NDIM);

  for(int i=0; i < 3; ++i)
    {
    size[i] = this->BoxLength/static_cast<float>(decompSize[i]);
    float m_corner_grid_alive = pos[i]*(this->NDIM/decompSize[i]);
    min[i] = m_corner_grid_alive*m_grid2phys_pos;
    } // END for
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::ClearParticlePositions(float **pos)
{
  free(pos[0]);
  free(pos);
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::PackageParticlePositions(
          SimulationParticles *particles, float ***positions)
{
  if( particles->NumParticles == 0 )
    {
    return;
    }

  *positions = (float **)malloc(sizeof(float*));
  (*positions)[0] = (float *)malloc(3*particles->NumParticles*sizeof(float));

  for( int i=0; i < particles->NumParticles; ++i )
    {
    (*positions)[0][i*3]   = particles->X[i];
    (*positions)[0][i*3+1] = particles->Y[i];
    (*positions)[0][i*3+2] = particles->Z[i];
    } // END for all particles
}

} /* namespace cosmotk */
