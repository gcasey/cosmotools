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



TessVoidFinderAnalysisTool::TessVoidFinderAnalysisTool()
{
  this->Name           = "TESS";
  this->Communicator   = MPI_COMM_NULL;
  this->MinVol         = .0001;
  this->MaxVol         = -1.0;
  this->GhostFactor    = 5.0;
  this->CellSize       = 1.0;
  this->Initialized    = false;
  this->TimeStatistics = new double[MAX_TIMES];
}

//------------------------------------------------------------------------------
TessVoidFinderAnalysisTool::~TessVoidFinderAnalysisTool()
{
  tess_finalize();
  if( this->TimeStatistics != NULL )
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

  int numBlocksPerProcess = 1;
  tess_init_diy_initialized(
      numBlocksPerProcess,this->CellSize,this->GhostFactor,
      this->MinVol,this->MaxVol,this->Communicator,this->TimeStatistics);

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
