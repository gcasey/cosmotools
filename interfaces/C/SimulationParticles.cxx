#include "SimulationParticles.h"

#include <cstdio> // For NULL

namespace cosmologytools
{

SimulationParticles::SimulationParticles()
{
  this->X  = this->Y  = this->Z  = NULL;
  this->VX = this->VY = this->VZ = NULL;
  this->GlobalIds           = NULL;
  this->HaloTags            = NULL;
  this->SubHaloTags         = NULL;
  this->NumParticles        = 0;
  this->TimeStep            = -1;
  this->RedShift            = 0.0;
  this->AreHalosComputed    = false;
  this->AreSubHalosComputed = false;
}

//-----------------------------------------------------------------------------
SimulationParticles::~SimulationParticles()
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

//-----------------------------------------------------------------------------
void SimulationParticles::AllocateHaloAndSubHaloArrays()
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

  this->AreHalosComputed    = false;
  this->AreSubHalosComputed = false;
}

} /* namespace cosmogolytools */
