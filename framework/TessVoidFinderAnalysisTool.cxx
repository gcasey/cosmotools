#include "TessVoidFinderAnalysisTool.h"

#include "diy.h"
#include "tess.h"

#include <cassert>

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
  this->MinVol = this->MaxVol = -1.0;
  this->GhostFactor = 5.0;
}

//------------------------------------------------------------------------------
TessVoidFinderAnalysisTool::~TessVoidFinderAnalysisTool()
{
  // TODO: implement this
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
  float **positions = NULL;
  this->PackageParticlePositions(particles, positions);
  // STEP 2: initialize tess, i.e., call tess_init()
  // TODO: implement this

  // STEP 3: execute tess on the given particle dataset, i.e., call tess()
  // TODO: implement this


  // STEP 4: finalize tess
  this->ClearParticlePositions(positions);
  tess_finalize();
}

//------------------------------------------------------------------------------
void TessVoidFinderAnalysisTool::WriteOutput()
{
  if( !this->GenerateOutput )
    {
    return;
    }

  // TODO: write output from tess. If tess is finalized, how to we do I/O here?
}

//------------------------------------------------------------------------------
std::string TessVoidFinderAnalysisTool::GetInformation()
{
  return(this->GetBasicInformation());
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
