/**
 * @class Halo
 * @brief A light-weight object to store the information for a single halo and
 * also provide functionality for performing operations between two halos,
 * e.g., intersection.
 */
#ifndef HALO_H_
#define HALO_H_

#include "CosmologyToolsMacros.h"

#include <iostream> // For ostream
#include <set> // For STL set

namespace cosmotk
{


class Halo
{
public:

  /**
   * @brief Default constructor
   */
  Halo();

  /**
   * @brief Custom constructor
   * @param Tag the tag of the halo
   * @param TimeStep the time-step the halo is
   * @param redShift the corresponding red-shift
   * @param cntr the halo-center
   * @param vel the velocity of the halo
   * @param particleIds the particle IDs of the halo
   * @param N the total number of particles
   */
  Halo(int Tag, int TimeStep, REAL redShift,
        POSVEL_T cntr[3], POSVEL_T vel[3],
        ID_T *particleIds, int N);

  /**
   * @brief Destructor.
   */
  virtual ~Halo();

  /**
   * @brief Intersects this halo instance with another halo
   * @param h the halo to intersect with
   * @return N the number of halos that intersect between this halo
   * instance and the halo pointed to by h.
   * @note iff N==0 this halo does not intersect with another halo.
   */
  int Intersect(Halo *h);

  /**
   * @brief Gets the hash code of this halo instance.
   * @return h a hashcode for this halo instance.
   */
  std::string GetHashCode();

  /**
   * @brief Prints a halo to the given C++ output stream
   * @param os the output stream object
   * @note Used for debugging.
   */
  void Print(std::ostream &os);

  int Tag;                      // The tag/ID of the halo
  int TimeStep;                 // The time-step of this halo
  REAL Redshift;                // The corresponding red-shift of the halo
  POSVEL_T Center[3];           // The halo-center
  POSVEL_T AverageVelocity[3];  // The average velocity of the halo
  std::set< ID_T > ParticleIds; // The global particle IDs of the halo
};

} /* namespace cosmotk */
#endif /* HALO_H_ */
