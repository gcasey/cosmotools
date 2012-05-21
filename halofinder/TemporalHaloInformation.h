/**
 * @class TemporalHaloInformation
 * @brief An object to hold two instances of HaloDataInformation at two different
 * time-steps. Besides acting as a storage facility for halo information, the
 * class also provide functionality for properly updating the temporal halo
 * information by internally management the current and previous pointers. The
 * current pointer will always point to the most recent halo information that
 * was computed while the previous pointer holds the results from k-timesteps
 * ago, where k is the user-supplied frequency.
 *
 * @see HaloInformation ForwardHaloTracker
 */
#ifndef TEMPORALHALOINFORMATION_H_
#define TEMPORALHALOINFORMATION_H_

#include "CosmologyToolsMacros.h"
#include <cstdlib> // For NULL

namespace cosmologytools
{

// Forward declarations
class HaloDataInformation;

class TemporalHaloInformation
{
public:
  TemporalHaloInformation();
  virtual ~TemporalHaloInformation();

  // Inline Get & Set macros
  GetNSetMacro(Current,HaloDataInformation*);
  GetNSetMacro(Previous,HaloDataInformation*);

  /**
   * @brief Updates this temporal halo information instance with halo data
   * at a new timestep. This method ensures proper execution and management of
   * Current and Previous halo information internally.
   * @param haloInformation pointer to the haloInformation object.
   */
  void Update(HaloDataInformation *haloInformation);

  /**
   * @brief Checks if the temporal information is complete. The temporal info
   * is complete if we have both Current and Previous halo information.
   * @return True iff (this->GetCurrent()!=NULL && this->GetPrevious()!=NULL)
   * holds. Otherwise, false.
   */
  bool IsComplete() {
    return( (this->Current!=NULL) && (this->Previous!=NULL)); };

protected:
  HaloDataInformation *Current;
  HaloDataInformation *Previous;

private:
  DISABLE_COPY_AND_ASSIGNMENT(TemporalHaloInformation);
};

} /* namespace cosmogolytools */
#endif /* TEMPORALHALOINFORMATION_H_ */
