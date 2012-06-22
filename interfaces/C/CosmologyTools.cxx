#include "CosmologyTools.h"

#include "SimulationParticles.h"

#include <iostream>
#include <cassert>

class SimulationParticles;

//------------------------------------------------------------------------------
void CosmologyInit(MPI_Comm *comm)
{
  CosmoToolsManager = new CosmologyToolsManager();
  if( comm != NULL )
    {
    CosmoToolsManager->SetCommunicator( *comm );
    }
  CosmoToolsManager->Barrier();
}

//------------------------------------------------------------------------------
void CosmologyFinit(MPI_Fint *fcomm)
{
  CosmoToolsManager = new CosmologyToolsManager();
  if( fcomm != NULL )
    {
    MPI_Comm comm     = MPI_Comm_f2c(*fcomm);
    CosmoToolsManager->SetCommunicator( comm );
    }
  CosmoToolsManager->Barrier();
}

//------------------------------------------------------------------------------
void CosmologyEnableVis()
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  CosmoToolsManager->SetEnableVis( true );
}

//------------------------------------------------------------------------------
void CosmologyDisableVis()
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  CosmoToolsManager->SetEnableVis( false );
}

//------------------------------------------------------------------------------
void CosmologySetParticles(
    INTEGER *tstep, REAL *redshift,
    REAL *px, REAL *py, REAL *pz,
    REAL *vx, REAL *vy, REAL *vz,
    INTEGER *GlobalParticlesIds,
    INTEGER *NumberOfParticles)
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  CosmoToolsManager->SetParticles(
      *tstep, *redshift,
      px, py,pz, vx,vy,vz,
      GlobalParticlesIds,
      *NumberOfParticles );
}

//------------------------------------------------------------------------------
void CosmologyGetHaloIds(INTEGER *haloTags)
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  assert("pre: user-supplied haloTags buffer is NULL" && (haloTags != NULL) );

  cosmologytools::SimulationParticles *myParticles =
      CosmoToolsManager->GetParticles();
  assert("pre: particle data-structure is NULL!" && (myParticles != NULL) );
  assert("pre: particles haloTags are NULL!" &&
          (myParticles->HaloTags != NULL) );

  for( int p=0; p < myParticles->NumParticles; ++p )
    {
    haloTags[ p ] = myParticles->HaloTags[ p ];
    } // END for all particles
}

//------------------------------------------------------------------------------
void CosmologyGetSubHaloIds(INTEGER *subHaloTags)
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  assert("pre: user-supplied subHaloTags is NULL" && (subHaloTags != NULL) );

 cosmologytools::SimulationParticles *myParticles =
     CosmoToolsManager->GetParticles();
 assert("pre: particle data-structure is NULL!" && (myParticles != NULL) );
 assert("pre: particles haloTags are NULL!" &&
         (myParticles->SubHaloTags != NULL) );

 for( int p=0; p < myParticles->NumParticles; ++p )
   {
   subHaloTags[ p ] = myParticles->SubHaloTags[ p ];
   } // END for all particles
}

//------------------------------------------------------------------------------
void CosmologySetMemoryLayout(int* memoryLayout)
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  assert("pre: undefined memory layout" &&
         (*memoryLayout >= 0) &&
         (*memoryLayout < MemoryLayout::NUMBER_OF_LAYOUTS) );

  CosmoToolsManager->SetLayout( *memoryLayout );
}

//------------------------------------------------------------------------------
void CosmologySetTrackerFrequency(int *frequency)
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  CosmoToolsManager->SetHaloTrackingFrequency( *frequency );
}

//------------------------------------------------------------------------------
void CosmologySetHaloFinder( int *haloFinder )
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  CosmoToolsManager->SetHaloFinder( *haloFinder );
}

//------------------------------------------------------------------------------
void CosmologyTrackHalos()
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  CosmoToolsManager->TrackHalos();
}

//------------------------------------------------------------------------------
void CosmologyFindHalos()
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  CosmoToolsManager->FindHalos();
}

//------------------------------------------------------------------------------
void CosmologyFinalize()
{
  assert("pre: CosmoToolsManager is NULL" && (CosmoToolsManager != NULL) );
  MPI_Comm comm = CosmoToolsManager->GetCommunicator();
  delete CosmoToolsManager;
  MPI_Barrier( comm );
}