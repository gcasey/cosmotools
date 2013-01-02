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

#include "diy.h" // For DIY_Datatype

/* TODO: for now, we explicitely define a max number of halos */
#define MAX_PARTICLES_IN_HALO 5

struct DIYHaloItem {
  int Tag;
  int TimeStep;
  REAL Redshift;
  POSVEL_T Center[3];
  POSVEL_T AverageVelocity[3];
  ID_T HaloParticles[MAX_PARTICLES_IN_HALO];
  int NumberOfHaloParticles;
};


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
   * @brief Construct a halo instance from the given DIYHaloItem
   * @param halo an instance of DIYHaloItem
   * @pre halo != NULL
   */
  Halo( DIYHaloItem *halo);

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

  /**
   * @brief Get a pointer to the corresponding DIY halo data type.
   * @return ptr Pointer to a registered DIY datatype.
   * @post ptr != NULL
   */
  DIY_Datatype* GetDIYHaloType();

  /**
   * @brief Populates a corresponding instance of DIYHaloItem
   * @param halo the DIYHaloItem instance to populate
   * @pre halo != NULL
   */
  void GetDIYHaloItem(DIYHaloItem *halo);

  /**
   * @brief Registers a DIY data-type to represent a halo object.
   * @note User is responsible for calling DIY_Destroy_datatype to
   * properly release the registered object.
   */
 static void CreateDIYHaloType(DIY_Datatype *dtype);


  int Tag;                      // The tag/ID of the halo
  int TimeStep;                 // The time-step of this halo
  REAL Redshift;                // The corresponding red-shift of the halo
  POSVEL_T Center[3];           // The halo-center
  POSVEL_T AverageVelocity[3];  // The average velocity of the halo
  std::set< ID_T > ParticleIds; // The global particle IDs of the halo


private:
  DIY_Datatype* DIYHaloType;

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
  void InitHalo(
          int Tag, int TimeStep, REAL redShift,
          POSVEL_T cntr[3], POSVEL_T vel[3],
          ID_T *particleIds, int N);
};

} /* namespace cosmotk */
#endif /* HALO_H_ */
