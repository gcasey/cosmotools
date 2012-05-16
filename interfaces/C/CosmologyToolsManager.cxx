#include "CosmologyToolsManager.h"

#include "HaloFinders.h"
#include "SimulationParticles.h"

#include <cassert>

namespace cosmologytools
{

CosmologyToolsManager::CosmologyToolsManager()
{
  this->Layout                = MemoryLayout::ROWMAJOR;
  this->Communicator          = MPI_COMM_WORLD;
  this->EnableVis             = false;
  this->HaloTrackingFrequency = 5;
  this->Particles             = new SimulationParticles();
  this->HaloFinder            = HaloFinders::COSMO;
}

//-----------------------------------------------------------------------------
CosmologyToolsManager::~CosmologyToolsManager()
{
  if( this->Particles != NULL )
    {
    delete this->Particles;
    }
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::SetParticles(
    INTEGER tstep, REAL redshift,
    REAL *px, REAL *py, REAL *pz,
    REAL *vx, REAL *vy, REAL *vz,
    INTEGER *GlobalParticlesIds,
    INTEGER NumberOfParticles)
{
  assert("pre: internal particles data is NULL" && (this->Particles!=NULL));

  this->Particles->TimeStep     = tstep;
  this->Particles->RedShift     = redshift;
  this->Particles->X            = px;
  this->Particles->Y            = py;
  this->Particles->Z            = pz;
  this->Particles->VX           = vx;
  this->Particles->VY           = vy;
  this->Particles->VZ           = vz;
  this->Particles->GlobalIds    = GlobalParticlesIds;
  this->Particles->NumParticles = NumberOfParticles;

  this->Particles->AllocateHaloAndSubHaloArrays();
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::TrackHalos()
{
 // TODO: implement this
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::FindHalos()
{
  // TODO: implement this
}

} /* namespace cosmologytools */
