#include "CosmologyToolsManager.h"

#include "HaloFinders.h"
#include "SimulationParticles.h"
#include "ForwardHaloTracker.h"

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
  this->HaloTracker           = new ForwardHaloTracker();
}

//-----------------------------------------------------------------------------
CosmologyToolsManager::~CosmologyToolsManager()
{
  if( this->Particles != NULL )
    {
    delete this->Particles;
    }

  if( this->HaloTracker != NULL )
    {
    delete this->HaloTracker;
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

  this->HaloTracker->RegisterParticles(
      tstep,redshift,px,py,pz,
      vx,vy,vz,NULL,NULL,
      GlobalParticlesIds,NULL,NULL,
      NumberOfParticles);
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::TrackHalos()
{
  this->HaloTracker->UpdateMergerTree(this->Particles->TimeStep);
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::FindHalos()
{
  // TODO: implement this
}

} /* namespace cosmologytools */
