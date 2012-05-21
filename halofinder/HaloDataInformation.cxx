#include "HaloDataInformation.h"

namespace cosmologytools
{

HaloDataInformation::HaloDataInformation()
{
  this->RedShift = 0.0;
  this->TimeStep = 0;
  this->GlobalIds.resize( 0 );
  this->HaloTags.resize( 0 );
}

//-----------------------------------------------------------------------------
HaloDataInformation::~HaloDataInformation()
{
  this->GlobalIds.clear();
  this->HaloTags.clear();
}

} /* namespace cosmogolytools */
