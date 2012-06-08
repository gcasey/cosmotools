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
  this->Communicator         = MPI_COMM_WORLD;
  this->Frequency            = 5;
  this->RL                   = 100;
  this->Overlap              = 5;
  this->LinkingLength        = 0.2;
  this->UseExplicitTimeSteps = false;
  this->TemporalHaloData     = new TemporalHaloInformation();
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
void ForwardHaloTracker::SetExplicitTrackerTimeSteps(
    INTEGER *tsteps, const INTEGER N)
{
  this->UseExplicitTimeSteps = true;
  for( int i=0; i < N; ++i )
    {
    this->TimeSteps.insert( tsteps[i] );
    }
}

//------------------------------------------------------------------------------
bool ForwardHaloTracker::IsTrackerTimeStep(const int tstep)
{
  if( this->UseExplicitTimeSteps )
    {
    if( this->TimeSteps.find( tstep ) != this->TimeSteps.end() )
      {
      return true;
      }
    else
      {
      return false;
      }
    }
  else
    {
    if( (tstep % this->Frequency) == 0)
      {
      return true;
      }
    else
      {
      return false;
      }
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
  if( this->IsTrackerTimeStep( tstep ) )
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
  if( this->IsTrackerTimeStep(tstep) && this->TemporalHaloData->IsComplete() )
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

  // STEP 0: Get halo-finder properties
  int numberOfHalos = hfinder->getNumberOfHalos();
  int *fofHalos     = hfinder->getHalos();
  int *fofHaloCount = hfinder->getHaloCount();
  int *fofHaloList  = hfinder->getHaloList();

  // STEP 1: Construct FOFHaloProperties object
  cosmologytools::FOFHaloProperties* fof =
      new cosmologytools::FOFHaloProperties();
  fof->setHalos(numberOfHalos,fofHalos,fofHaloCount,fofHaloList);
  fof->setParameters("",this->RL,this->Overlap,this->LinkingLength);

  // STEP 1: Filter out the halos
  std::vector< int > extractedHalos;
  for(int halo=0; halo < numberOfHalos; ++halo )
    {
    if( fofHaloCount[halo] >= this->PMin )
      {
      extractedHalos.push_back( halo );
      } // END if the halo is within the pmin threshold
    } // END for all halos

  // STEP 2: Extract particle halo ids and particle ids
  for(int halo=0; halo < static_cast<int>(extractedHalos.size()); ++halo )
    {
    int internalHaloIdx = extractedHalos[ halo ];
    } // END for all extracted halos

}

} /* namespace cosmogolytools */
