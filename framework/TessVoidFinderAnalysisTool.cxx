#include "TessVoidFinderAnalysisTool.h"

#include "diy.h"
#include "tess.h"
#include "voronoi.h"

#include <cassert>
#include <sstream>

#include <mpi.h>

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

namespace cosmotk
{

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
  this->PackageParticlePositions(particles, positions);

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

  // STEP 0: Ensure we are dealing with a cartesian communicator
  int topology = 0;
  MPI_Topo_test(this->Communicator,&topology);
  assert("pre: communicator is not cartesian!" && topology==MPI_CART);
  if( topology != MPI_CART )
    {
    std::cerr << "ERROR: tess is expecting a cartesian communicator!\n";
    MPI_Abort(this->Communicator,-1);
    }

  // STEP 1: Compute number of blocks per process and total number of blocks
  // Currently, numBlocksPerProcess = 1 and totalNumberOfBlocks is given by
  // numBlocksPerProcess * numRanks.
  int numBlocksPerProcess = 1;
  int totalBlocks = 0;
  int numRanks = 0;
  int rank = 0;
  MPI_Comm_rank(this->Communicator,&rank);
  MPI_Comm_size(this->Communicator,&numRanks);
  totalBlocks = numBlocksPerProcess * numRanks;

  // STEP 2: We are dealing with a 3-D domain and with structured topology, so,
  // each block has 26 neighbors.
  const int DIMENSION        = 3;
  const int NUM_OF_NEIGHBORS = 26;

  // STEP 3: Get the Block bounds
  bb_t bb;
  float min[DIMENSION], size[DIMENSION];

// TODO: I need to get this information somehow
// Domain::rL_local_alive(size);
// Domain::corner_phys_alive(min);

  for( int i=0; i < DIMENSION; ++i )
    {
    bb.min[i] = min[i];
    bb.max[i] = min[i] + size[i];
    }

  // STEP 4: Get decomposition size (number of blocks in each dimension)
  int decomp_size[DIMENSION];
// TODO: Need a way to get the decomposition size, from the communicator perhaps?
// Partition::getDecompSize(decomp_size);

  // data overall extents
  // assume all blocks are same size (as mine)
  // assume 0,0,0 is the overall data minimum corner
  float data_mins[DIMENSION], data_maxs[DIMENSION];
  for (int i = 0; i < DIMENSION; i++) {
    data_mins[i] = 0.0;
    data_maxs[i] = data_mins[i] + size[i] * decomp_size[i];
  }

  // STEP 5: Get neighbors
  int neigh_gids[NUM_OF_NEIGHBORS];
// TODO: Need a way to get this information
//  Partition::getNeighbors(neigh_gids);
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

  // STEP 7: Get wrap
  int wrap = (this->Periodic)? 1:0;

  // STEP 8: call tess_init
  tess_init(numBlocksPerProcess,totalBlocks,gids,&bb,neighbors,num_neighbors,
            this->CellSize,this->GhostFactor,data_mins,data_maxs,wrap,
            this->MinVol,this->MaxVol,this->Communicator,
            this->TimeStatistics);

  // STEP 9: Return all dynamically allocated memory
  delete[] neighbors[0];
  delete[] neighbors;

  // STEP 10: Set initialized to true
  this->Initialized = true;
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::ClearParticlePositions(float **pos)
{
  free(pos[0]);
  free(pos);
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::PackageParticlePositions(
          SimulationParticles *particles, float **positions)
{
  if( particles->NumParticles == 0 )
    {
    return;
    }

  positions = (float **)malloc(sizeof(float*));
  positions[0] = (float *)malloc(3*particles->NumParticles*sizeof(float));

  for( int i=0; i < particles->NumParticles; ++i )
    {
    positions[0][i*3]   = particles->X[i];
    positions[0][i*3+1] = particles->Y[i];
    positions[0][i*3+2] = particles->Z[i];
    } // END for all particles
}

} /* namespace cosmotk */
