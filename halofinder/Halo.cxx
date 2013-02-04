#include "Halo.h"

#include <algorithm>
#include <cassert>
#include <sstream>
#include <vector>
#include <cstdio>
#include <cstddef>

namespace cosmotk
{

//-----------------------------------------------------------------------------
void Halo::CreateDIYHaloType(DIY_Datatype *dtype)
{
  struct map_block_t halo_map[5] = {
   {DIY_INT,     OFST, 1, offsetof(struct DIYHaloItem, Tag)},
   {DIY_INT,     OFST, 1, offsetof(struct DIYHaloItem, TimeStep)},
   {DIY_REAL_T,   OFST, 1, offsetof(struct DIYHaloItem, Redshift)},
   {DIY_POSVEL_T, OFST, 3, offsetof(struct DIYHaloItem, Center)},
   {DIY_POSVEL_T, OFST, 3, offsetof(struct DIYHaloItem, AverageVelocity)},
  };
  DIY_Create_struct_datatype(0, 5, halo_map, dtype);
}

//-----------------------------------------------------------------------------
void Halo::CreateDIYHaloParticleType(DIY_Datatype *dtype)
{
  struct map_block_t halo_part_map[3] = {
   {DIY_INT,  OFST, 1, offsetof(struct DIYHaloParticleItem, Tag)},
   {DIY_INT,  OFST, 1, offsetof(struct DIYHaloParticleItem, TimeStep)},
   {DIY_ID_T,  OFST, 1, offsetof(struct DIYHaloParticleItem, HaloParticleID)},
  };
  DIY_Create_struct_datatype(0, 3, halo_part_map, dtype);
}

//-----------------------------------------------------------------------------
Halo::Halo()
{
  this->Tag = -1;
  this->TimeStep = 0;
  this->Redshift = 0.;
  this->HaloMass = 0.;
  this->Center[0] = this->Center[1] = this->Center[2] = 0.;
  this->AverageVelocity[0] =
  this->AverageVelocity[1] =
  this->AverageVelocity[2] = 0.;
}

//-----------------------------------------------------------------------------
Halo::Halo( DIYHaloItem *halo )
{
  assert("pre: DIYHaloItem is NULL!" && (halo != NULL) );
  this->InitHalo(
        halo->Tag,
        halo->TimeStep,
        halo->Redshift,
        halo->Center,
        halo->AverageVelocity,
        NULL,
        0 );
}

//-----------------------------------------------------------------------------
Halo::Halo(
    int Tag, int TimeStep, REAL redShift,
    POSVEL_T cntr[3], POSVEL_T vel[3],
    ID_T *particleIds, int N)
{
  assert("pre: the particle IDs array should not be NULL" &&
          (particleIds != NULL) );
  this->InitHalo(Tag,TimeStep,redShift,cntr,vel,particleIds,N);
}

//-----------------------------------------------------------------------------
Halo::~Halo()
{
  this->ParticleIds.clear();
}

//-----------------------------------------------------------------------------
void Halo::InitHalo(
    int Tag, int TimeStep, REAL redShift,
    POSVEL_T cntr[3], POSVEL_T vel[3],
    ID_T *particleIds, int N)
{
  this->Tag      = Tag;
  this->TimeStep = TimeStep;
  this->Redshift = redShift;
  for( int i=0; i < 3; ++i )
    {
    this->Center[i]          = cntr[i];
    this->AverageVelocity[i] = vel[i];
    }

  if(particleIds==NULL)
    {
    return;
    }

  for( int i=0; i < N; ++i )
    {
    this->ParticleIds.insert( particleIds[i] );
    }
}

//-----------------------------------------------------------------------------
void Halo::SetHaloParticles(ID_T *particleIds, int N)
{
  this->ParticleIds.clear();
  if( particleIds == NULL )
    {
    return;
    }

  for(int i=0; i < N; ++i)
    {
    this->ParticleIds.insert(particleIds[i]);
    } // END for
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

//-----------------------------------------------------------------------------
std::string Halo::GetHashCode()
{
  return( Halo::GetHashCodeForHalo(this->Tag,this->TimeStep));
}

//-----------------------------------------------------------------------------
std::string Halo::GetHashCodeForHalo(int tag, int timestep)
{
  std::ostringstream oss;
  oss << timestep << "." << tag;
  return( oss.str() );
}

//-----------------------------------------------------------------------------
void Halo::Print(std::ostream &os)
{
  os << std::endl;
  os << "================================\n";
  os << "HALO TAG: "    << this->Tag      << std::endl;
  os << "TSTEP: "       << this->TimeStep << std::endl;
  os << "RED-SHIFT: "   << this->Redshift << std::endl;
  os << "HALO Center: ";
  for(int i=0; i < 3; os << this->Center[i++] << " ");
  os << std::endl;
  os << "HALO Average Velocity: ";
  for(int i=0; i < 3; os << this->AverageVelocity[i++] << " ");
  os << std::endl;

  os << "HALO Particles: " << std::endl << "\t";
  std::set< ID_T >::iterator iter = this->ParticleIds.begin();
  for( ; iter != this->ParticleIds.end(); ++iter )
    {
    os << *iter << " ";
    } // END for all halo particles
  os << std::endl;
}

//-----------------------------------------------------------------------------
void Halo::GetDIYHaloItem(DIYHaloItem *halo)
{
  assert("pre: DIYHaloItem instance is NULL" && (halo != NULL) );

  halo->Tag      = this->Tag;
  halo->TimeStep = this->TimeStep;
  halo->Redshift = this->Redshift;

  for( int i=0; i < 3; ++i )
    {
    halo->Center[i]          = this->Center[i];
    halo->AverageVelocity[i] = this->AverageVelocity[i];
    }
}

//-----------------------------------------------------------------------------
void Halo::GetDIYHaloParticleItemsVector(
          std::vector<DIYHaloParticleItem> &haloParticles)
{
  haloParticles.resize(this->ParticleIds.size());
  std::set< ID_T >::iterator iter = this->ParticleIds.begin();
  for(int idx=0; iter != this->ParticleIds.end(); ++iter, ++idx)
    {
    haloParticles[ idx ].Tag            = this->Tag;
    haloParticles[ idx ].TimeStep       = this->TimeStep;
    haloParticles[ idx ].HaloParticleID = *iter;
    } // END for all particles
}

} /* namespace cosmotk */
