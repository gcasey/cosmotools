#include "ParallelHaloMergerTree.h"

// cosmotools includes
#include "Halo.h"
#include "HaloNeighborExchange.h"
#include "HaloType.h"
#include "MPIUtilities.h"

// DIY communication sub-strate
#include "diy.h"

#include <fstream>
#include <sstream>

namespace cosmotk
{

ParallelHaloMergerTree::ParallelHaloMergerTree()
{
  this->Communicator     = MPI_COMM_NULL;
  this->CurrentIdx       = this->PreviousIdx = -1;
  this->NeighborExchange = new HaloNeighborExchange();
}

//------------------------------------------------------------------------------
ParallelHaloMergerTree::~ParallelHaloMergerTree()
{
  this->TemporalHalos.clear();

  if( this->NeighborExchange != NULL )
    {
    delete this->NeighborExchange;
    }
}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::UpdateMergerTree(
    const int t1, Halo *haloSet1, const int M,
    const int t2, Halo *haloSet2, const int N,
    DistributedHaloEvolutionTree *t)
{
  assert("pre: t1 < t2" && (t1 < t2) );
  assert("pre: halo evolution tree is NULL!" && (t != NULL) );
  assert("pre: HaloNeighborExchange object is NULL!" &&
          (this->NeighborExchange != NULL) );

  this->NeighborExchange->SetCommunicator( this->Communicator );

  // STEP 0: Exchange halos -- This step set's up the TemporalHalos data-struct
  if( this->TemporalHalos.size() == 0 )
    {
    // This is the first time that we run the algorithm, so we must communicate
    // both time-steps
    this->TemporalHalos.resize(2);
    this->PreviousIdx = 0;
    this->CurrentIdx  = 1;

    this->NeighborExchange->ExchangeHalos(
        haloSet1,M,this->TemporalHalos[this->PreviousIdx]);

    this->NeighborExchange->ExchangeHalos(
        haloSet2,N,this->TemporalHalos[this->CurrentIdx]);
    }
  else
    {
    assert("pre: TemoralHalos.size()==2" && (this->TemporalHalos.size()==2));

    // Instead of copying vectors, we just swap indices to the TemporalHalos
    std::swap(this->CurrentIdx,this->PreviousIdx);

    // Exchange halos at the current time-step. Note!!!! halos at the
    // previous time-step were exchanged at the previous time-step.

    this->NeighborExchange->ExchangeHalos(
        haloSet2,N,this->TemporalHalos[this->CurrentIdx]);
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

  // STEP 4: Print diagnostics for this interval
  if( this->Verbose )
    {
    MPIUtilities::Printf(
     this->Communicator,"=============================\n");
    MPIUtilities::Printf(
     this->Communicator,"TIMESTEP INTERVAL [%d,%d]\n",t1,t2);
    MPIUtilities::Printf(
     this->Communicator,"TOTAL NUMBER OF BIRTHS: %d\n",
     this->GetTotalNumberOfBirths());
    MPIUtilities::Printf(
     this->Communicator,"TOTAL NUMBER OF RE-BIRTHS: %d\n",
     this->GetTotalNumberOfRebirths());
    MPIUtilities::Printf(
     this->Communicator,"TOTAL NUMBER OF MERGERS: %d\n",
     this->GetTotalNumberOfMerges());
    MPIUtilities::Printf(
     this->Communicator,"TOTAL NUMBER OF SPLITS: %d\n",
     this->GetTotalNumberOfSplits());
    MPIUtilities::Printf(
     this->Communicator,"TOTAL NUMBER OF DEATHS: %d\n",
     this->GetTotalNumberOfDeaths());
    MPIUtilities::Printf(
     this->Communicator,"=============================\n");
    }

  // STEP 5: Handle Death events
  this->HandleDeathEvents( t );

  // STEP 6: Barrier synchronization
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
      HaloType::SetType(zombie->HaloTypeMask,HaloType::ZOMBIE);
      zombie->Tag = (-1)*zombie->Tag;
      }

    // More sanity checks
    assert("pre: Node not flagged as a zombie!" &&
           (HaloType::IsType(zombie->HaloTypeMask,HaloType::ZOMBIE)));
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

} /* namespace cosmotk */
