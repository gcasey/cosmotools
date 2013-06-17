/**
 * @class HaloDataInformation
 * @brief A light-weight class to hold the halo-information at a given time-step
 * i.e., the time-step and corresponding red-shift, the global particles IDs
 * and the corresponding halo tags.
 *
 * @see TemporalHaloInformation ForwardHaloTracker
 */
#ifndef HALODATAINFORMATION_H_
#define HALODATAINFORMATION_H_

#include "CosmoToolsMacros.h"
#include <vector> // For STL vector

#include "Halo.h"

namespace cosmologytools
{

class HaloDataInformation
{
public:
  HaloDataInformation();
  virtual ~HaloDataInformation();

  REAL RedShift;
  int TimeStep;
  int NumberOfHalos;

  std::vector<cosmotk::Halo> Halos;

private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloDataInformation);

};

} /* namespace cosmogolytools */
#endif /* HALODATAINFORMATION_H_ */
