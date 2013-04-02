#include "ParallelHaloMergerTree.h"

// cosmotools includes
#include "Halo.h"
#include "MPIUtilities.h"

// DIY communication sub-strate
#include "diy.h"

#include <fstream>
#include <sstream>

namespace cosmotk
{

ParallelHaloMergerTree::ParallelHaloMergerTree()
{
  this->CurrentIdx = this->PreviousIdx = -1;
}

//------------------------------------------------------------------------------
ParallelHaloMergerTree::~ParallelHaloMergerTree()
{
  this->TemporalHalos.clear();
}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::UpdateMergerTree(
    const int t1, Halo *haloSet1, const int M,
    const int t2, Halo *haloSet2, const int N,
    DistributedHaloEvolutionTree *t)
{
  assert("pre: t1 < t2" && (t1 < t2) );
  assert("pre: halo evolution tree is NULL!" && (t != NULL) );


  // STEP 0: Exchange halos -- This step set's up the TemporalHalos data-struct
  if( this->TemporalHalos.size() == 0 )
    {
    // This is the first time that we run the algorithm, so we must communicate
    // both time-steps
    this->TemporalHalos.resize(2);
    this->PreviousIdx = 0;
    this->CurrentIdx  = 1;

    cosmotk::MPIUtilities::Printf(
        this->Communicator,"Exchanging cached halos...");
    this->ExchangeHalos(haloSet1,M,this->TemporalHalos[this->PreviousIdx]);
    cosmotk::MPIUtilities::Printf(
        this->Communicator,"[DONE]\n");

    cosmotk::MPIUtilities::Printf(
            this->Communicator,"Exchanging halos at current timestep...");
    this->ExchangeHalos(haloSet2,N,this->TemporalHalos[this->CurrentIdx]);
    cosmotk::MPIUtilities::Printf(
            this->Communicator,"[DONE]\n");
    }
  else
    {
    assert("pre: TemoralHalos.size()==2" && (this->TemporalHalos.size()==2));

    // Instead of copying vectors, we just swap indices to the TemporalHalos
    std::swap(this->CurrentIdx,this->PreviousIdx);

    // Exchange halos at the current time-step. Note!!!! halos at the
    // previous time-step were exchanged at the previous time-step.
    this->ExchangeHalos(haloSet2,N,this->TemporalHalos[this->CurrentIdx]);
    }

  // STEP 1: Register halos
  this->RegisterHalos(
      t1,&this->TemporalHalos[this->PreviousIdx][0],
          this->TemporalHalos[this->PreviousIdx].size(),
      t2,&this->TemporalHalos[this->CurrentIdx][0],
          this->TemporalHalos[this->CurrentIdx].size());

  // STEP 2: Compute merger-tree
  this->ComputeMergerTree();
  //this->PrintMatrix();

  // STEP 3: Update the halo-evolution tree
  this->UpdateHaloEvolutionTree( t );

  // STEP 4: Handle Death events
  this->HandleDeathEvents( t );

  // STEP 5: Barrier synchronization
  this->Barrier();
}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::HandleDeathEvents(
        DistributedHaloEvolutionTree *t)
{
  assert("pre: halo evolution tree is NULL!" && (t != NULL) );

  std::set< int >::iterator iter = this->DeadHalos.begin();
  for(;iter != this->DeadHalos.end(); ++iter)
    {
    int haloIdx = *iter;

    // Sanity checks
    assert("pre: Dead haloIdx is out-of-bounds" &&
        (haloIdx >= 0) &&
        (haloIdx < this->TemporalHalos[this->PreviousIdx].size()) );

    unsigned char bitmask;
    MergerTreeEvent::Reset(bitmask);
    MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::DEATH);

    // Copy the halo
    this->TemporalHalos[this->CurrentIdx].push_back(
          this->TemporalHalos[this->PreviousIdx][haloIdx]);
    int zombieIdx = this->TemporalHalos[this->CurrentIdx].size()-1;

    // Get pointers to the zombie halo at the current timestep and the
    // source halo at the previous timestep, i.e., the halo that died.
    Halo *zombie     = &this->TemporalHalos[this->CurrentIdx][zombieIdx];
    Halo *sourceHalo = &this->TemporalHalos[this->PreviousIdx][haloIdx];

    // Get a pointer to a reference halo at the current timestep so that we
    // can update the timestep and red-shift information.
    Halo *refHalo    = &this->TemporalHalos[this->CurrentIdx][0];
    zombie->TimeStep = refHalo->TimeStep;
    zombie->Redshift = refHalo->Redshift;
    if( zombie->Count == 0 )
      {
      // this halo just died for the first time
      zombie->HaloType = ZOMBIEHALO;
      zombie->Tag = (-1)*zombie->Tag;
      }

    // More sanity checks
    assert("pre: Node not flagged as a zombie!" &&
           (zombie->HaloType==ZOMBIEHALO));
    assert("pre: zombies must have a negative tag!" && (zombie->Tag < 0));
    zombie->Count++;

    // Insert the zombie halo in to the tree
    this->InsertHalo(zombie,bitmask,t);

    // Create a link between the source halo and zombie halo
    t->LinkHalos(sourceHalo->GetHashCode(),zombie->GetHashCode());
    } // END for all dead halos
}

//------------------------------------------------------------------------------
int ParallelHaloMergerTree::GetTotalNumberOfBirths()
{
  int total = 0;
  MPI_Allreduce(
      &this->NumberOfBirths,&total,1,MPI_INT,MPI_SUM,this->Communicator);
  return( total );
}

//------------------------------------------------------------------------------
int ParallelHaloMergerTree::GetTotalNumberOfRebirths()
{
  int total = 0;
  MPI_Allreduce(
      &this->NumberOfRebirths,&total,1,MPI_INT,MPI_SUM,this->Communicator);
  return( total );
}

//------------------------------------------------------------------------------
int ParallelHaloMergerTree::GetTotalNumberOfMerges()
{
  int total = 0;
  int local = this->MergeHalos.size();
  MPI_Allreduce(&local,&total,1,MPI_INT,MPI_SUM,this->Communicator);
  return( total );
}

//------------------------------------------------------------------------------
int ParallelHaloMergerTree::GetTotalNumberOfSplits()
{
  int total = 0;
  int local = this->SplitHalos.size();
  MPI_Allreduce(&local,&total,1,MPI_INT,MPI_SUM,this->Communicator);
  return( total );
}

//------------------------------------------------------------------------------
int ParallelHaloMergerTree::GetTotalNumberOfDeaths()
{
  int total = 0;
  int local = this->DeadHalos.size();
  MPI_Allreduce(&local,&total,1,MPI_INT,MPI_SUM,this->Communicator);
  return( total );
}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::ExchangeHalos(
      cosmotk::Halo *halos, const int N,
      std::vector<cosmotk::Halo>& globalHalos)
{
  assert("pre: N >= 0" && (N >= 0) );

  HaloHashMap haloHash;
  this->ExchangeHaloInfo( halos, N, haloHash );
  this->ExchangeHaloParticles( halos, N, haloHash );

  globalHalos.resize( N+haloHash.size() );

  int haloIdx = 0;
  for(; haloIdx < N; ++haloIdx )
    {
    globalHalos[haloIdx] = halos[haloIdx];
    } // END for all local halos

  HaloHashMap::iterator iter = haloHash.begin();
  for(;iter != haloHash.end(); ++iter, ++haloIdx )
    {
    globalHalos[haloIdx] = iter->second;
    } // END for all neighboring halos
}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::ExchangeHaloInfo(
        cosmotk::Halo *halos, const int N, HaloHashMap& haloHash)
{
  assert("pre: N >= 0" && (N >= 0) );
  assert("pre: haloHash.empty()" && haloHash.empty() );

  // STEP 0: Enqueue halo information
  HaloInfo haloInfo;
  for( int hidx=0; hidx < N; ++hidx )
    {
    halos[ hidx ].GetHaloInfo( &haloInfo );
    DIY_Enqueue_item_all(
        0, 0, (void *)&haloInfo, NULL, sizeof(HaloInfo), NULL);
    } // END for all halos

  // STEP 1: Neighbor exchange
  int nblocks          = 1;
  void ***rcvHalos     = new void**[nblocks];
  int *numHalosReceived = new int[nblocks];
  DIY_Exchange_neighbors(
      0,rcvHalos,numHalosReceived,1.0,&cosmotk::Halo::CreateDIYHaloInfoType);

  // STEP 2: Unpack received halos and store them in the halo hash by
  // a halo hash code.
  HaloInfo *rcvHaloItem = NULL;
  for(int i=0; i < nblocks; ++i)
    {
    for( int j=0; j < numHalosReceived[i]; ++j )
      {
      rcvHaloItem = (struct HaloInfo *)rcvHalos[i][j];
      cosmotk::Halo h(rcvHaloItem);
      haloHash[ h.GetHashCode() ] = h;
      } // END for all received halos of this block
    } // END for all blocks

  // STEP 3: Clean up
  DIY_Flush_neighbors(
      0,rcvHalos,numHalosReceived,&cosmotk::Halo::CreateDIYHaloInfoType);
  delete [] numHalosReceived;

}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::ExchangeHaloParticles(
      cosmotk::Halo *halos, const int N, HaloHashMap& haloHash)
{
  // STEP 0: For each halo, enqueue its halo particle IDs
  std::vector<HaloParticle> haloParticles;
  for( int hidx=0; hidx < N; ++hidx )
    {
    halos[hidx].GetHaloParticlesVector(haloParticles);
    for( int pidx=0; pidx < haloParticles.size(); ++pidx )
      {
      DIY_Enqueue_item_all(
          0,0,
          (void*)&haloParticles[pidx],
          NULL,
          sizeof(HaloParticle),
          NULL);
      } // END for all halo particles
    } // END for all halos

  // STEP 1: Neighbor exchange
  int nblocks = 1;
  void ***rcvHalos = new void**[nblocks];
  int *numHalosReceived = new int[nblocks];
  DIY_Exchange_neighbors(
      0,rcvHalos,numHalosReceived,1.0,
      &cosmotk::Halo::CreateDIYHaloParticleType);

  // STEP 2: Unpack data to neighbor halos in the haloHash
  HaloParticle *haloParticle = NULL;
  for( int i=0; i < nblocks; ++i )
    {
    for( int j=0; j < numHalosReceived[i]; ++j )
      {
      haloParticle = (struct HaloParticle*)rcvHalos[i][j];
      std::string hashCode =
          cosmotk::Halo::GetHashCodeForHalo(
              haloParticle->Tag,haloParticle->TimeStep);
      assert(haloHash.find(hashCode)!=haloHash.end());
      haloHash[hashCode].ParticleIds.insert(haloParticle->HaloParticleID);
      } // END for all received halos of this block
    } // END for all blocks

  // STEP 3: Clean up
  DIY_Flush_neighbors(
    0,rcvHalos,numHalosReceived,&cosmotk::Halo::CreateDIYHaloParticleType);
  delete [] numHalosReceived;
}

} /* namespace cosmotk */
