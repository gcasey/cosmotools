#include "CosmologyToolsManager.h"

#include "HaloFinders.h"
#include "SimulationParticles.h"
#include "ForwardHaloTracker.h"

#include <cstdlib>
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

  this->TimeSteps = NULL;
  this->NumTimeSteps = 0;
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
void CosmologyToolsManager::SetAnalysisTimeSteps(
        INTEGER *timeSteps, INTEGER numTimeSteps)
{
  assert("pre: timesteps is NULL!" && (timeSteps != NULL) );
  assert("pre: numTimeSteps > 0" && (numTimeSteps > 0) );

  this->TimeSteps = timeSteps;
  this->NumTimeSteps = numTimeSteps;

  // Set explicit time-steps for the halo-tracker
  this->HaloTracker->SetExplicitTrackerTimeSteps(
      this->TimeSteps,this->NumTimeSteps);

  // TODO: Set explicit time-steps also for the halo-finder ?
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::SetParticles(
    INTEGER tstep, REAL redshift,
    POSVEL_T *px, POSVEL_T *py, POSVEL_T *pz,
    POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz,
    REAL *mass, POTENTIAL_T *potential,
    ID_T *GlobalParticlesIds,
    MASK_T *mask,
    STATUS_T *state,
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
  this->Particles->Mass         = mass;
  this->Particles->Potential    = potential;
  this->Particles->GlobalIds    = GlobalParticlesIds;
  this->Particles->Mask         = mask;
  this->Particles->State        = state;
  this->Particles->NumParticles = NumberOfParticles;
  this->Particles->AllocateHaloAndSubHaloArrays();

  this->HaloTracker->RegisterParticles(
      tstep,redshift,px,py,pz,
      vx,vy,vz,mass,potential,
      GlobalParticlesIds,mask,state,
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
