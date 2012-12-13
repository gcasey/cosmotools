#include "Halo.h"

#include <algorithm>
#include <cassert>
#include <vector>

namespace cosmotk
{

//-----------------------------------------------------------------------------
Halo::Halo()
{
  this->Tag = -1;
  this->TimeStep = 0;
  this->Redshift = 0.;
  this->Center[0] = this->Center[1] = this->Center[2] = 0.;
  this->AverageVelocity[0] =
  this->AverageVelocity[1] =
  this->AverageVelocity[2] = 0.;
}

//-----------------------------------------------------------------------------
Halo::Halo(
    int Tag, int TimeStep, REAL redShift,
    POSVEL_T cntr[3], POSVEL_T vel[3],
    ID_T *particleIds, int N)
{
  assert("pre: the particle IDs array should not be NULL" &&
          (particleIds != NULL) );

  this->Tag      = Tag;
  this->TimeStep = TimeStep;
  this->Redshift = redShift;
  for( int i=0; i < 3; ++i )
    {
    this->Center[i]          = cntr[i];
    this->AverageVelocity[i] = vel[i];
    }
  for( int i=0; i < N; ++i )
    {
    this->ParticleIds.insert( particleIds[i] );
    }

}

//-----------------------------------------------------------------------------
Halo::~Halo()
{
  this->ParticleIds.clear();
}

//-----------------------------------------------------------------------------
int Halo::Intersect(Halo *h)
{
  assert("pre: halo to intersect with should not be NULL" && (h != NULL) );

  int overlap = 0;

  std::vector<ID_T> overlapParticles;
  std::set_intersection(
      this->ParticleIds.begin(),this->ParticleIds.end(),
      h->ParticleIds.begin(),h->ParticleIds.end(),
      std::inserter(overlapParticles, overlapParticles.end()));

  overlap = overlapParticles.size();
  overlapParticles.clear();
  return( overlap );
}

} /* namespace cosmotk */
