#include "HaloDataInformation.h"

namespace cosmologytools
{

HaloDataInformation::HaloDataInformation()
{
  this->RedShift = 0.0;
  this->TimeStep = 0;
}

//-----------------------------------------------------------------------------
HaloDataInformation::~HaloDataInformation()
{
  this->Halos.clear();
}

} /* namespace cosmogolytools */
