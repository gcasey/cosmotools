#include "CosmologyToolsManager.h"
#include "HaloFinders.h"

#include <cassert>

namespace cosmologytools
{

//------------------------------------------------------------------------------
//  SIMULATION PARTICLES CLASS
//------------------------------------------------------------------------------

/**
 * @class Particles
 * @brief An object to hold a reference to the simulation arrays for the
 * particle data. Currently the particle data consists of a 3-component
 */
class SimulationParticles
{
public:

  /**
   * @brief Default constructor
   */
  SimulationParticles()
    {
    this->X  = this->Y  = this->Z  = NULL;
    this->VX = this->VY = this->VZ = NULL;
    this->GlobalIds    = NULL;
    this->HaloTags     = NULL;
    this->SubHaloTags  = NULL;
    this->NumParticles = 0;
    this->TimeStep     = -1;
    this->RedShift     = 0.0;
    }

  /**
   * @brief Destructor
   */
  ~SimulationParticles()
    {
    // NOTE: these are pointers to the simulation's memory space, thus,
    // the co-processing library will not delete these pointers.
    this->X  = this->Y  = this->Z  = NULL;
    this->VX = this->VY = this->VZ = NULL;
    this->GlobalIds    = NULL;
    this->NumParticles = 0;

    if( this->HaloTags != NULL )
      {
      delete [] this->HaloTags;
      }

    if( this->SubHaloTags != NULL )
      {
      delete [] this->SubHaloTags;
      }
    }

  /**
   * @brief Allcates the halo and subhalo arrays
   */
  void AllocateHaloAndSubHaloArrays()
  {
    if( this->HaloTags != NULL )
      {
      delete [] this->HaloTags;
      }
    this->HaloTags = new INTEGER[ this->NumParticles ];

    if( this->SubHaloTags != NULL )
      {
      delete [] this->SubHaloTags;
      }
    this->SubHaloTags = new INTEGER[ this->NumParticles ];
  }

  INTEGER TimeStep;     /** the current timestep */
  REAL RedShift;        /** the red-shift */
  REAL *X;              /** x-component of the particles position vector  */
  REAL *Y;              /** y-component of the particles position vector  */
  REAL *Z;              /** z-component of the particles position vector  */
  REAL *VX;             /** vx-component of the particles position vector */
  REAL *VY;             /** vy-component of the particles position vector */
  REAL *VZ;             /** vz-component of the particles position vector */
  INTEGER *GlobalIds;   /** list of global ids */
  INTEGER NumParticles; /** the total number of particles */

  INTEGER *HaloTags;    /** List of halo ids for each particle (computed) */
  INTEGER *SubHaloTags; /** List of subhalo ids for each particles (computed) */

private:
  DISABLE_COPY_AND_ASSIGNMENT(SimulationParticles);
};

//------------------------------------------------------------------------------
//  CosmologyToolsManager implementation
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------

CosmologyToolsManager::CosmologyToolsManager()
{
  this->Layout                = MemoryLayout::ROWMAJOR;
  this->Communicator          = MPI_COMM_WORLD;
  this->EnableVis             = false;
  this->HaloTrackingFrequency = 5;
  this->Particles             = new SimulationParticles();
  this->HaloFinder            = HaloFinders::COSMO;
}

//-----------------------------------------------------------------------------
CosmologyToolsManager::~CosmologyToolsManager()
{
  if( this->Particles != NULL )
    {
    delete this->Particles;
    }
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::SetParticles(
    INTEGER tstep, REAL redshift,
    REAL *px, REAL *py, REAL *pz,
    REAL *vx, REAL *vy, REAL *vz,
    INTEGER *GlobalParticlesIds,
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
  this->Particles->GlobalIds    = GlobalParticlesIds;
  this->Particles->NumParticles = NumberOfParticles;

  this->Particles->AllocateHaloAndSubHaloArrays();
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::TrackHalos()
{
 // TODO: implement this
}

//-----------------------------------------------------------------------------
void CosmologyToolsManager::FindHalos()
{
  // TODO: implement this
}

} /* namespace cosmologytools */
