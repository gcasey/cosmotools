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
  this->NumberOfParticles  = 0;
}

//------------------------------------------------------------------------------
ForwardHaloTracker::~ForwardHaloTracker()
{
  if (this->TemporalHaloData != NULL)
  {
  delete this->TemporalHaloData;
  }
  this->NumberOfParticles = 0;
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
    REAL* vx, REAL* vy, REAL *vz,
    REAL* mass, REAL* potential,
    INTEGER* id, INTEGER* mask, INTEGER* state,
    INTEGER N)
{
  // Sanity checks
  assert("pre: x particles positions array is NULL!" && (px != NULL) );
  assert("pre: y particles positions array is NULL!" && (py != NULL));
  assert("pre: z particles positions array is NULL!" && (pz != NULL) );
  assert("pre: vx paritlces velocities array is NULL!" && (vx != NULL) );
  assert("pre: vy particles velocities array is NULL!" && (vy != NULL) );
  assert("pre: vz particles velocities array is NULL!" && (vz != NULL) );
  assert("pre: mass particles array is NULL!" && (mass != NULL) );
  assert("pre: mask particles array is NULL!" && (mask != NULL) );
  assert("pre: state particles arrray is NULL!" && (state != NULL) );

  if( this->IsTrackerTimeStep( tstep ) )
    {
    this->NumberOfParticles = N;

    CosmoHaloFinderP *haloFinder = new CosmoHaloFinderP;

    // Create vectors from the pointers since the halo-finder accepts vectors
    // NOTE: these are slow...halo-finder ought to work with pointers, my 2c.
    this->Px.assign(px, px+N);
    this->Py.assign(py, py+N);
    this->Pz.assign(pz, pz+N);
    this->Vx.assign(vx, vx+N);
    this->Vy.assign(vy, vy+N);
    this->Vz.assign(vz, vz+N);
    this->Potential.assign(potential,potential+N);
    this->Mass.assign(mass,mass+N);
    this->Id.assign(id,id+N);
    this->Mask.assign(mask,mask+N);
    this->State.assign(state,state+N);

    haloFinder->setParticles(
        &(this->Px),&(this->Py),&(this->Pz),
        &(this->Vx),&(this->Vy),&(this->Vz),
        &(this->Potential),
        &(this->Id),
        &(this->Mask),
        &(this->State)
        );

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
  fof->setParticles(
      &(this->Px),&(this->Py),&(this->Pz),
      &(this->Vx),&(this->Vy),&(this->Vz),
      &(this->Mass),
      &(this->Potential),
      &(this->Id),
      &(this->Mask),
      &(this->State)
      );

  // STEP 2: Filter out the halos within the PMIN threshold
  std::vector< int > extractedHalos;
  for(int halo=0; halo < numberOfHalos; ++halo )
    {
    if( fofHaloCount[halo] >= this->PMin )
      {
      extractedHalos.push_back( halo );
      } // END if the halo is within the pmin threshold
    } // END for all halos


  hinfo->Allocate( extractedHalos.size() );

  // STEP 3: Extract particle halo ids and particle ids
  for(int halo=0; halo < static_cast<int>(extractedHalos.size()); ++halo )
    {
    int internalHaloIdx     = extractedHalos[ halo ];
    int NumParticlesInHalo  = fofHaloCount[ internalHaloIdx ];

    // TODO: we should change the way we access the halo-finder information
    // here. We need to have a more efficient and intuitive API.
    REAL *xlocHalo = new REAL[ NumParticlesInHalo ];
    REAL *ylocHalo = new REAL[ NumParticlesInHalo ];
    REAL *zlocHalo = new REAL[ NumParticlesInHalo ];
    REAL *xVelHalo = new REAL[ NumParticlesInHalo ];
    REAL *yVelHalo = new REAL[ NumParticlesInHalo ];
    REAL *zVelHalo = new REAL[ NumParticlesInHalo ];
    REAL *massHalo = new REAL[ NumParticlesInHalo ];
    int  *id       = new int[ NumParticlesInHalo ];
    int *actualIdx = new int[ NumParticlesInHalo ];

    fof->extractInformation(
        internalHaloIdx,actualIdx,
        xlocHalo,ylocHalo,zlocHalo,
        xVelHalo,yVelHalo,zVelHalo,
        massHalo,id);

    for(int hpidx=0; hpidx < NumParticlesInHalo; ++hpidx )
      {
      hinfo->GlobalIds.push_back( id[ actualIdx[hpidx] ] );
      hinfo->HaloTags.push_back( internalHaloIdx );
      } // END for all particles in halo

    delete [] xlocHalo;
    delete [] ylocHalo;
    delete [] zlocHalo;
    delete [] xVelHalo;
    delete [] yVelHalo;
    delete [] zVelHalo;
    delete [] massHalo;
    delete [] id;
    delete [] actualIdx;
    } // END for all extracted halos


  // STEP 4: Delete FOFHaloProperties
  delete fof;
}

} /* namespace cosmogolytools */
