#include "ForwardHaloTracker.h"

#include "CosmoHaloFinderP.h"
#include "FOFHaloProperties.h"
#include "HaloDataInformation.h"
#include "HaloFinders.h"
#include "Partition.h"
#include "TemporalHaloInformation.h"

#include <cassert>
#include <vector>

namespace cosmologytools {

ForwardHaloTracker::ForwardHaloTracker()
{
  this->BoxLength         = 0.0;
  this->NG                = 5;
  this->NDIM              = 0;
  this->PMIN              = 200;
  this->LinkingLength     = 0.2;
  this->NumberOfParticles = 0;
  this->Communicator      = MPI_COMM_NULL;

  this->TemporalHaloData  = new TemporalHaloInformation();
  this->HaloEvolutionTree = new cosmotk::DistributedHaloEvolutionTree();
  this->HaloMergerTree    = new cosmotk::ParallelHaloMergerTree();
}

//------------------------------------------------------------------------------
ForwardHaloTracker::~ForwardHaloTracker()
{
  if( this->TemporalHaloData != NULL )
    {
    delete this->TemporalHaloData;
    }

  if( this->HaloEvolutionTree != NULL )
    {
    this->HaloEvolutionTree->WriteTree();
    delete this->HaloEvolutionTree;
    }

  if( this->HaloMergerTree != NULL )
    {
    delete this->HaloMergerTree;
    }
  this->NumberOfParticles = 0;
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::RegisterParticles(
    const INTEGER tstep, const REAL redShift,
    POSVEL_T* px, POSVEL_T* py, POSVEL_T*pz,
    POSVEL_T* vx, POSVEL_T* vy, POSVEL_T*vz,
    POSVEL_T* mass, POTENTIAL_T* potential,
    ID_T* id, MASK_T* mask, STATUS_T* state,
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

  this->TimeStep = tstep;
  this->RedShift = redShift;

  this->Px = px; this->Py = py; this->Pz = pz;
  this->Vx = vx; this->Vy = vy; this->Vz = vz;
  this->Mass              = mass;
  this->Potential         = potential;
  this->Mask              = mask;
  this->State             = state;
  this->NumberOfParticles = N;
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::TrackHalos()
{
  // STEP 0: Execute the halo-finder
  cosmologytools::Partition::initialize(this->Communicator);
  cosmologytools::CosmoHaloFinderP *haloFinder=
     new cosmologytools::CosmoHaloFinderP();
  this->ExecuteHaloFinder(haloFinder);

  // STEP 1: Get halos at this time-step. Note, the pointer to currentHaloData
  // is managed by the TemporalHaloData object.
  HaloDataInformation *currentHaloData = new HaloDataInformation();
  this->GetHaloInformation(currentHaloData,haloFinder);

  // STEP 2: Update temporal halo information
  this->TemporalHaloData->Update(currentHaloData);

  // STEP 3: Update merger-tree
  this->UpdateMergerTree();
  this->Barrier();
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::ExecuteHaloFinder(
        CosmoHaloFinderP* haloFinder)
{
  assert("pre: HaloFinder is NULL" && (haloFinder != NULL) );

  // STEP 0: Set the parameters for the halo-finder
  haloFinder->setParameters(
      "",this->BoxLength,this->NG,this->NDIM,this->PMIN,
      this->LinkingLength);

  // STEP 1: Set the particle input dataset
  haloFinder->setParticles(
      this->Px,this->Py,this->Pz,
      this->Vx,this->Vy,this->Vz,
      this->Potential,this->Id,
      this->Mask,this->State,
      this->NumberOfParticles);

  // STEP 2: Run the halo-finder in parallel
  haloFinder->executeHaloFinder();

  // STEP 3: Merge results across ranks
  haloFinder->collectHalos(false /*clearSerial*/);
  haloFinder->mergeHalos();
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::UpdateMergerTree()
{
  if( this->TemporalHaloData->IsComplete() )
    {
    HaloDataInformation* current  = this->TemporalHaloData->GetCurrent();
    HaloDataInformation* previous = this->TemporalHaloData->GetPrevious();
    this->HaloMergerTree->UpdateMergerTree(
        current->TimeStep,&current->Halos[0],current->NumberOfHalos,
        previous->TimeStep,&previous->Halos[0],previous->NumberOfHalos,
        this->HaloEvolutionTree);
    } // END if update merger-tree
}

//------------------------------------------------------------------------------
void ForwardHaloTracker::GetHaloInformation(
    HaloDataInformation* haloData, CosmoHaloFinderP* haloFinder)
{
  assert("pre: halo information object is NULL" && haloData != NULL);
  assert("pre: hfinder information object is NULL" && haloFinder != NULL);

  // STEP 0: Get halo-finder properties
  int numberOfHalos = haloFinder->getNumberOfHalos();
  int *fofHalos     = haloFinder->getHalos();
  int *fofHaloCount = haloFinder->getHaloCount();
  int *fofHaloList  = haloFinder->getHaloList();

  // STEP 1: Construct FOFHaloProperties object
  cosmologytools::FOFHaloProperties* fof =
      new cosmologytools::FOFHaloProperties();
  fof->setHalos(numberOfHalos,fofHalos,fofHaloCount,fofHaloList);
  fof->setParameters("",this->BoxLength,this->NG,this->LinkingLength);
  fof->setParticles(
      this->NumberOfParticles,
      this->Px,this->Py,this->Pz,
      this->Vx,this->Vy,this->Vz,
      this->Mass,
      this->Id );

  // STEP 2: Extract halo information
  haloData->NumberOfHalos = 0;
  haloData->TimeStep      = this->TimeStep;
  haloData->RedShift      = this->RedShift;

  cosmotk::Halo myHalo;

  for( int halo=0; halo < numberOfHalos; ++halo)
    {
    if( fofHaloCount[halo] >= this->PMIN )
      {
      myHalo.Tag = haloData->NumberOfHalos;
      haloData->NumberOfHalos++;

      myHalo.TimeStep = this->TimeStep;
      myHalo.Redshift = this->RedShift;

      ID_T* haloParticles = new ID_T[fofHaloCount[halo]];
      fof->extractHaloParticleIds(halo,haloParticles);
      for( int hpidx=0; hpidx < fofHaloCount[halo]; ++hpidx)
        {
        myHalo.ParticleIds.insert( haloParticles[hpidx] );
        } // END for all halo particles
      delete [] haloParticles;

      fof->FOFHaloPosition(halo,myHalo.Center);
      fof->FOFHaloVelocity(halo,myHalo.AverageVelocity);

      haloData->Halos.push_back(myHalo);
      } // END if halo is big enough
    } // END for all halos

  // STEP 3: Delete FOFHaloProperties
  delete fof;
}

} /* namespace cosmogolytools */
