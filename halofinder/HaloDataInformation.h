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

#include "CosmologyToolsMacros.h"
#include <vector> // For STL vector

namespace cosmologytools
{

class HaloDataInformation
{
public:
  HaloDataInformation();
  virtual ~HaloDataInformation();

  /**
   * @brief Allocates internal data-structures
   * @param N the number of halos
   */
  void Allocate(int N)
    {
    this->NumberOfHalos = N;
    this->HaloCenter.resize( 3*N, 0.0 );
    this->HaloAverageVelocity.resize( 3*N, 0.0 );
    this->HaloVelocityDispersion.resize( N, 0.0 );
    this->HaloMass.resize( N, 0.0 );
    }

  REAL RedShift;
  int TimeStep;
  int NumberOfHalos;

  std::vector< REAL > HaloCenter;
  std::vector< REAL > HaloVelocityDispersion;
  std::vector< REAL > HaloAverageVelocity;
  std::vector< REAL > HaloMass;
  std::vector< INTEGER > GlobalIds;
  std::vector< INTEGER > HaloTags;

private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloDataInformation);

};

} /* namespace cosmogolytools */
#endif /* HALODATAINFORMATION_H_ */
