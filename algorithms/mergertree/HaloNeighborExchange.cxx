#include "HaloNeighborExchange.h"

// CosmologyTools includes
#include "Halo.h"
#include "HaloType.h"
#include "MPIUtilities.h"

// DIY
#include "diy.h"


namespace cosmotk
{

HaloNeighborExchange::HaloNeighborExchange()
{
  this->Communicator = MPI_COMM_NULL;

}

//------------------------------------------------------------------------------
HaloNeighborExchange::~HaloNeighborExchange()
{
  this->EnqueuedHalos.clear();
}

//------------------------------------------------------------------------------
void HaloNeighborExchange::EnqueueHalo(Halo* h)
{
  assert("pre: input halo is NULL!" && (h != NULL) );
  this->EnqueuedHalos.push_back( *h );

}

//------------------------------------------------------------------------------
void HaloNeighborExchange::ExchangeHalos(
    std::vector<Halo>& exchangedHalos,bool inclusive)
{
  this->ExchangeHalos(
      &this->EnqueuedHalos[0],
      this->EnqueuedHalos.size(),exchangedHalos,
      inclusive);

  this->EnqueuedHalos.clear();
}

//------------------------------------------------------------------------------
void HaloNeighborExchange::ExchangeHalos(
    Halo *localHalos, const int N,
    std::vector<Halo>& exchangedHalos, bool inclusive)
{
  // STEP 0: Get halos from neighbors. This is done in two steps. First, the
  // halo information of each halo is exchanged. This step also initializes the
  // 'neighborHalos' data-structure. Second, the particle IDs of each halo
  // are exchanged.
  HaloHashMap neighborHalos;
  this->ExchangeHaloInformation(localHalos,N,neighborHalos);
  this->ExchangeHaloParticles(localHalos,N,neighborHalos);

  // STEP 1: Allocate output vector size to store both local and neighboring
  // halos.
  int size = (inclusive)? N+neighborHalos.size() : neighborHalos.size();
  exchangedHalos.resize(size);

  // STEP 2: Put local halos in the output vector iff inclusive
  int haloIdx = 0;

  if( inclusive )
    {
    for(; haloIdx < N; ++haloIdx )
      {
      exchangedHalos[ haloIdx ]= localHalos[ haloIdx ];
      } // END for all local halos
    } // END if inclusive

  // STEP 3: Put neighboring halos in the output vector
  HaloHashMap::iterator iter = neighborHalos.begin();
  for(; iter != neighborHalos.end(); ++iter, ++haloIdx )
    {
    exchangedHalos[ haloIdx ] = iter->second;
    assert("pre: neighbor halo should be marked as a ghost!" &&
        (HaloType::IsType(exchangedHalos[haloIdx].HaloTypeMask,
                          HaloType::GHOST)) );
    } // END for all neighboring halos

  // STEP 4: Clean up
  neighborHalos.clear();
}

//------------------------------------------------------------------------------
void HaloNeighborExchange::ExchangeHaloInformation(
      Halo *localHalos, const int N, HaloHashMap& neighborHalos )
{
  assert("pre: neighborHalos.empty()" && neighborHalos.empty() );

  // STEP 0: Enqueue HaloInfo objects to send to neighbors
  HaloInfo hinfo;
  for( int hidx=0; hidx < N; ++hidx )
    {
    localHalos[ hidx ].GetHaloInfo( &hinfo );
    DIY_Enqueue_item_all(
       0, 0, (void*)&hinfo, NULL, sizeof(HaloInfo), NULL);
    } // END for all local halos

  // STEP 1: Allocate receive buffer and exchange data with neighbors
  int nblocks           = 1;
  void ***rcvHaloInfo   = new void **[nblocks];
  int *numHalosReceived = new int [nblocks];
  DIY_Exchange_neighbors(
     0,rcvHaloInfo,numHalosReceived,1.0,&Halo::CreateDIYHaloInfoType);

  // STEP 2: Unpack neighbor halo information
  HaloInfo *hinforcv = NULL;
  for( int i=0; i < nblocks; ++i )
    {
    for( int j=0; j < numHalosReceived[i]; ++j )
      {
      hinforcv = (struct HaloInfo*)rcvHaloInfo[i][j];
      Halo nei(hinforcv);

      assert("pre: duplicate neighbor halo detected!" &&
          neighborHalos.find(nei.GetHashCode()) == neighborHalos.end() );

      neighborHalos[ nei.GetHashCode() ] = nei;
      } // END for all received halos of this block
    } // END for all blocks

  // STEP 3: Clean up
  DIY_Flush_neighbors(
     0,rcvHaloInfo,numHalosReceived,&Halo::CreateDIYHaloInfoType);
  delete [] numHalosReceived;

}

//------------------------------------------------------------------------------
void HaloNeighborExchange::ExchangeHaloParticles(
      Halo *localHalos, const int N, HaloHashMap& neighborHalos )
{
  // STEP 0: Enqueue HaloParticle objects to send to neighbors
  std::vector< HaloParticle > hparticles;
  for(int hidx=0; hidx < N; ++hidx)
    {
    localHalos[ hidx ].GetHaloParticlesVector( hparticles );
    for(int pidx=0; pidx < hparticles.size(); ++pidx )
      {
      DIY_Enqueue_item_all(
         0, 0, (void*)&hparticles[pidx], NULL, sizeof(HaloParticle),
         NULL);
      } // END for all halo particles
    } // END for all halos

  // STEP 1: Allocate receive buffer and exchange data with neighbors
  int nblocks              = 1;
  void ***rcvHaloParticles = new void**[nblocks];
  int *numParticlesRcvd    = new int[ nblocks];
  DIY_Exchange_neighbors(
    0,rcvHaloParticles,numParticlesRcvd,1.0,&Halo::CreateDIYHaloParticleType);

  // STEP 2:Unpack halo particles
  HaloParticle *haloParticle = NULL;
  for(int i=0; i < nblocks; ++i)
    {
    for(int j=0; j < numParticlesRcvd[i]; ++j)
      {
      haloParticle = (struct HaloParticle*)rcvHaloParticles[i][j];
      std::string hashCode =
          Halo::GetHashCodeForHalo(haloParticle->Tag,haloParticle->TimeStep);
      assert( neighborHalos.find(hashCode) != neighborHalos.end() );
      neighborHalos[hashCode].ParticleIds.insert(
          haloParticle->HaloParticleID);
      } // END for all particles of this block
    } // END for all blocks

  // STEP 3: Clean up
  DIY_Flush_neighbors(
    0, rcvHaloParticles, numParticlesRcvd, &Halo::CreateDIYHaloParticleType);
  delete [] numParticlesRcvd;
}

} /* namespace cosmotk */
