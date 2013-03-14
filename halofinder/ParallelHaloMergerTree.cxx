#include "ParallelHaloMergerTree.h"

#include "Halo.h"

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
void ParallelHaloMergerTree::PrintMatrix()
{
  int rank;
  MPI_Comm_rank(this->Communicator,&rank);

  std::ostringstream oss;
  oss << "HaloSimilarityMatrix_P" << rank
      << "_" << this->Timesteps[0] << "-" << this->Timesteps[1] << ".dat";

  std::ofstream ofs;
  ofs.open( oss.str().c_str() );
  assert( ofs.is_open() );

  int nrows = this->Sizes[0];
  int ncol  = this->Sizes[1];
  for(int row=0; row < nrows; ++row)
    {
    for(int col=0; col < ncol; ++col)
      {
      ofs << this->HaloSimilarityMatrix[row*ncol+col] << " ";
      } // END for all columns
    ofs << std::endl;
    } // END for all rows

  ofs.close();
}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::UpdateMergerTree(
    const int t1, Halo *haloSet1, const int M,
    const int t2, Halo *haloSet2, const int N,
    DistributedHaloEvolutionTree *t)
{
  assert("pre: t1 < t2" && (t1 < t2) );
  assert("pre: halos is NULL" && (haloSet1 != NULL) );
  assert("pre: halos is NULL" && (haloSet2 != NULL) );
  assert("pre: M > 0" && (M > 0) );
  assert("pre: N > 0" && (N > 0) );
  assert("pre: halo evolution tree is NULL!" && (t != NULL) );


  // STEP 0: Exchange halos -- This step set's up the TemporalHalos data-struct
  if( this->TemporalHalos.size() == 0 )
    {
    // This is the first time that we run the algorithm, so we must communicate
    // both time-steps
    this->TemporalHalos.resize(2);
    this->PreviousIdx = 0;
    this->CurrentIdx  = 1;
    this->ExchangeHalos(haloSet1,M,this->TemporalHalos[this->PreviousIdx]);
    this->ExchangeHalos(haloSet2,N,this->TemporalHalos[this->CurrentIdx]);
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
    assert("pre: Dead haloIdx is out-of-bounds" &&
        (haloIdx >= 0) &&
        (haloIdx < this->TemporalHalos[this->PreviousIdx].size()) );

    // Copy the halo
    this->TemporalHalos[this->CurrentIdx].push_back(
          this->TemporalHalos[this->PreviousIdx][haloIdx]);
    int zombieIdx = this->TemporalHalos[this->CurrentIdx].size()-1;

    // Update its attributes
    Halo *zombie     = &this->TemporalHalos[this->CurrentIdx][zombieIdx];
    Halo *sourceHalo = &this->TemporalHalos[this->PreviousIdx][haloIdx];
    Halo *refHalo    = &this->TemporalHalos[this->CurrentIdx][0];
    zombie->TimeStep = refHalo->TimeStep;
    zombie->Redshift = refHalo->Redshift;
    zombie->HaloType = ZOMBIEHALO;

    t->InsertNode( *zombie );
    t->CreateEdge(
         sourceHalo->GetHashCode(),
         zombie->GetHashCode(),
         100,
         MergerTreeEvent::DEATH);
    } // END for all dead halos
}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::ExchangeHalos(
      cosmotk::Halo *halos, const int N,
      std::vector<cosmotk::Halo>& globalHalos)
{
  assert("pre: halos is NULL!" && (halos != NULL) );
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
  assert("pre: halos is NULL!" && (halos != NULL) );
  assert("pre: N >= 0" && (N >= 0) );
  assert("pre: haloHash.empty()" && haloHash.empty() );

  // STEP 0: Enqueue halo information
  DIYHaloItem haloInfo;
  for( int hidx=0; hidx < N; ++hidx )
    {
    halos[ hidx ].GetDIYHaloItem( &haloInfo );
    DIY_Enqueue_item_all(
        0, 0, (void *)&haloInfo, NULL, sizeof(DIYHaloItem), NULL);
    } // END for all halos

  // STEP 1: Neighbor exchange
  int nblocks          = 1;
  void ***rcvHalos     = new void**[nblocks];
  int *numHalosReceived = new int[nblocks];
  DIY_Exchange_neighbors(
      0,rcvHalos,numHalosReceived,1.0,&cosmotk::Halo::CreateDIYHaloType);

  // STEP 2: Unpack received halos and store them in the halo hash by
  // a halo hash code.
  DIYHaloItem *rcvHaloItem = NULL;
  for(int i=0; i < nblocks; ++i)
    {
    for( int j=0; j < numHalosReceived[i]; ++j )
      {
      rcvHaloItem = (struct DIYHaloItem *)rcvHalos[i][j];
      cosmotk::Halo h(rcvHaloItem);
      haloHash[ h.GetHashCode() ] = h;
      } // END for all received halos of this block
    } // END for all blocks

  // STEP 3: Clean up
  DIY_Flush_neighbors(
      0,rcvHalos,numHalosReceived,&cosmotk::Halo::CreateDIYHaloType);
  delete [] numHalosReceived;

}

//------------------------------------------------------------------------------
void ParallelHaloMergerTree::ExchangeHaloParticles(
      cosmotk::Halo *halos, const int N, HaloHashMap& haloHash)
{
  // STEP 0: For each halo, enqueue its halo particle IDs
  std::vector<DIYHaloParticleItem> haloParticles;
  for( int hidx=0; hidx < N; ++hidx )
    {
    halos[hidx].GetDIYHaloParticleItemsVector(haloParticles);
    for( int pidx=0; pidx < haloParticles.size(); ++pidx )
      {
      DIY_Enqueue_item_all(
          0,0,
          (void*)&haloParticles[pidx],
          NULL,
          sizeof(DIYHaloParticleItem),
          NULL);
      } // END for all halo particles
    } // END for all halos

  // STEP 1: Neighbor exchange
  int nblocks = 1;
  void ***rcvHalos = new void**[nblocks];
  int *numHalosReceived = new int[nblocks];
  DIY_Exchange_neighbors(
      0,rcvHalos,numHalosReceived,1.0,&cosmotk::Halo::CreateDIYHaloParticleType);

  // STEP 2: Unpack data to neighbor halos in the haloHash
  DIYHaloParticleItem *haloParticle = NULL;
  for( int i=0; i < nblocks; ++i )
    {
    for( int j=0; j < numHalosReceived[i]; ++j )
      {
      haloParticle = (struct DIYHaloParticleItem*)rcvHalos[i][j];
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
