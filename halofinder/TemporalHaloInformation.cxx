#include "TemporalHaloInformation.h"
#include "HaloDataInformation.h"

#include <cstdlib>
#include <cassert>

namespace cosmologytools
{

//------------------------------------------------------------------------------
TemporalHaloInformation::TemporalHaloInformation()
{
  this->Current = this->Previous = NULL;
}

//------------------------------------------------------------------------------
TemporalHaloInformation::~TemporalHaloInformation()
{
  if( this->Current != NULL )
   {
   delete this->Current;
   this->Current = NULL;
   }

  if( this->Previous != NULL )
   {
   delete this->Previous;
   this->Previous = NULL;
   }
}

//------------------------------------------------------------------------------
void TemporalHaloInformation::Update(HaloDataInformation *haloInfo)
{
  assert("pre: haloInfo is NULL" && (haloInfo != NULL));

  if( this->Current != NULL )
   {
   if(this->Previous != NULL)
     {
     delete this->Previous;
     }
   this->Previous = this->Current;
   this->Current  = haloInfo;
   }
  else
   {
   // The first time this method is called, both Current and Previous are NULL.
   // We set the supplied data to the current pointer, and previous is set to
   // NULL here just for sanity.
   this->Current  = haloInfo;
   this->Previous = NULL;
   }
}

//------------------------------------------------------------------------------
bool TemporalHaloInformation::IsComplete()
{
  if( (this->Current != NULL)  &&
      (this->Previous != NULL) &&
      (this->Current->NumberOfHalos > 0) &&
      (this->Previous->NumberOfHalos > 0) )
      return true;
  return false;
}

} /* namespace cosmogolytools */
