#include "SimulationParticles.h"

#include <cstdio>  // For NULL
#include <cassert> // For assert

namespace cosmotk
{

SimulationParticles::SimulationParticles()
{
  this->X  = this->Y  = this->Z  = NULL;
  this->VX = this->VY = this->VZ = NULL;
  this->Mass      = NULL;
  this->Potential = NULL;
  this->GlobalIds = NULL;
  this->Mask      = NULL;
  this->State     = NULL;
  this->NumParticles = 0;
  this->TimeStep     = -1;
  this->RedShift     = 0.0;

}

//-----------------------------------------------------------------------------
SimulationParticles::~SimulationParticles()
{
  // NOTE: these are pointers to the simulation's memory space, thus,
  // the co-processing library will not delete these pointers.
  this->X  = this->Y  = this->Z  = NULL;
  this->VX = this->VY = this->VZ = NULL;
  this->Mass         = NULL;
  this->Potential    = NULL;
  this->GlobalIds    = NULL;
  this->Mask         = NULL;
  this->State        = NULL;
  this->NumParticles = 0;
}

//-----------------------------------------------------------------------------
void SimulationParticles::SetParticles(
    INTEGER tstep, REAL redShift,
    POSVEL_T *px, POSVEL_T *py, POSVEL_T *pz,
    POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz,
    POSVEL_T *mass, POTENTIAL_T *potential,
    ID_T *tags, MASK_T *mask, STATUS_T *status,
    ID_T N)
{
  assert("pre: NULL array supplied!" && (px != NULL) );
  assert("pre: NULL array supplied!" && (py != NULL) );
  assert("pre: NULL array supplied!" && (pz != NULL) );
  assert("pre: NULL array supplied!" && (vx != NULL) );
  assert("pre: NULL array supplied!" && (vy != NULL) );
  assert("pre: NULL array supplied!" && (vz != NULL) );
  assert("pre: NULL array supplied!" && (mass != NULL) );
  assert("pre: NULL array supplied!" && (potential != NULL) );
  assert("pre: NULL array supplied!" && (tags != NULL) );
  assert("pre: NULL array supplied!" && (mask != NULL) );
  assert("pre: NULL array supplied!" && (status != NULL) );

  this->TimeStep     = tstep;
  this->RedShift     = redShift;
  this->X            = px;
  this->Y            = py;
  this->Z            = pz;
  this->VX           = vx;
  this->VY           = vy;
  this->VZ           = vz;
  this->Mass         = mass;
  this->Potential    = potential;
  this->GlobalIds    = tags;
  this->Mask         = mask;
  this->State        = status;
  this->NumParticles = N;
}

} // END cosmotk namespace
