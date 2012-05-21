#include "ForwardHaloTracker.h"

#include "FOFHaloProperties.h"
#include "CosmoHaloFinderP.h"
#include "HaloFinders.h"
#include "HaloDataInformation.h"
#include "TemporalHaloInformation.h"

#include <cassert>
#include <vector>

namespace cosmologytools {

ForwardHaloTracker::ForwardHaloTracker()
{
  this->Communicator     = MPI_COMM_WORLD;
  this->Frequency        = 5;
  this->RL               = 100;
  this->Overlap          = 5;
  this->LinkingLength    = 0.2;
  this->TemporalHaloData = new TemporalHaloInformation();
}

//------------------------------------------------------------------------------
ForwardHaloTracker::~ForwardHaloTracker()
{
  if (this->TemporalHaloData != NULL)
  {
  delete this->TemporalHaloData;
  }
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::RegisterParticles(
    const int tstep, const double redShift,
    REAL* px, REAL* py, REAL *pz,
    REAL* vx, REAL* vy, REAL *vz, REAL* potential,
    INTEGER* id, INTEGER* mask, INTEGER* state,
    INTEGER N)
{
  if( (tstep % this->Frequency) == 0)
    {
    CosmoHaloFinderP *haloFinder = new CosmoHaloFinderP;

    // Create vectors from the pointers since the halo-finder accepts vectors
    std::vector< REAL > xLoc, yLoc, zLoc, xVel, yVel, zVel, potentialVector;
    std::vector< INTEGER > idVector, stateVector;
    std::vector< short unsigned int > maskVector;

    // NOTE: these are slow...halo-finder ought to work with pointers, my 2c.
    xLoc.assign(px, px+N);
    yLoc.assign(py, py+N);
    zLoc.assign(pz, pz+N);
    xVel.assign(vx, vx+N);
    yVel.assign(vy, vy+N);
    zVel.assign(vz, vz+N);
    potentialVector.assign(potential,potential+N);
    idVector.assign(id,id+N);
    maskVector.assign(mask,mask+N);
    stateVector.assign(state,state+N);

    haloFinder->setParticles(
        &xLoc,&yLoc,&zLoc,&xVel,&yVel,&zVel,
        &potentialVector,&idVector,&maskVector,&stateVector);

    this->Barrier();
    haloFinder->executeHaloFinder();
    this->Barrier();
    haloFinder->collectHalos();
    this->Barrier();
    haloFinder->mergeHalos();
    this->Barrier();

    // Extract the halo information data at this time-step. Note, these objects
    // are managed by the TemporalHaloInformation object. TemporalHaloData
    // de-allocates the HaloDataInformation once they go out-of-scope, that's
    // why the object is not de-allocated here.
    HaloDataInformation *hinfo = new HaloDataInformation();
    hinfo->RedShift = redShift;
    hinfo->TimeStep = tstep;
    this->GetHaloInformation(hinfo, haloFinder);

    // Update the Temporal halo-information
    this->TemporalHaloData->Update(hinfo);
    delete haloFinder;
    }

  // Synch all procs
  this->Barrier();
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::UpdateMergerTree(const int tstep)
{
  if( (tstep%this->Frequency==0) && this->TemporalHaloData->IsComplete() )
    {
    // TODO: update the merger trees here
    // Plug-in Jay's code here!
    }
  this->Barrier();
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::GetHaloInformation(
    HaloDataInformation* hinfo, CosmoHaloFinderP* hfinder)
{
  assert("pre: halo information object is NULL" && hinfo != NULL);
  assert("pre: hfinder information object is NULL" && hfinder != NULL);

  // TODO: Extract the global particle IDs and halo tags
}

} /* namespace cosmogolytools */
