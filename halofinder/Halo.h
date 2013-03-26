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
#include <set>      // For STL set
#include <vector>   // For STL vector

#include "diy.h" // For DIY_Datatype

/**
 * @struct HaloInfo
 * @brief Used to encapsulate static size halo information for communication
 * over DIY. This is also the same data-structure that is stored in the
 * merger-tree.
 */
struct HaloInfo {
  int Tag;
  int TimeStep;
  REAL Redshift;
  REAL HaloMass;
  POSVEL_T Center[3];
  POSVEL_T MeanCenter[3];
  POSVEL_T AverageVelocity[3];
  int DIYGlobalId;
};

/**
 * @struct DIYHaloParticleItem
 * @brief Used to encapsulate a *single* halo particle ID for communication
 * over DIY.
 */
struct DIYHaloParticleItem {
  int Tag;
  int TimeStep;
  ID_T HaloParticleID;
};

namespace cosmotk
{

enum HaloTypeEnum {
  NORMALHALO, // normal halos are halos that are owned by the given process
  GHOSTHALO,  // Ghost halos are halos that exist on other processes
  ZOMBIEHALO  // Zombie halos, are halos that existed in some previous time-step
              // and then seized to exist
};
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
   * @brief Construct a halo instance from the given HaloInfo
   * @param halo an instance of HaloInfo
   * @pre halo != NULL
   */
  Halo( HaloInfo *halo);

  /**
   * @brief Destructor.
   */
  virtual ~Halo();

  /**
   * @brief Sets the halo particles of this instance.
   * @param particleIds pointer to an array of particle IDs
   * @param N the number of particles.
   * @note particleIds==NULL iff N==0.
   */
  void SetHaloParticles(ID_T *particleIds, int N);

  /**
   * @brief Intersects this halo instance with another halo
   * @param h the halo to intersect with
   * @return N the percentage of particles of the given halo, h, that intersect
   * with this halo instance.
   * @note iff N==0 this halo does not intersect with another halo.
   */
  int Intersect(Halo *h);

  /**
   * @brief Returns the the number of particles in the halo
   * @return N the number of particles
   */
  int GetNumberOfParticles()
    {return static_cast<int>(this->ParticleIds.size());}

  /**
   * @brief Computes a hashcode for a halo with the given tag(ID) and timestep.
   * @param tag the tag or ID of the halo
   * @param timestep the timestep of the halo
   * @return h a hashcode for a halo with the given tag and timestep
   */
  static std::string GetHashCodeForHalo(int tag, int timestep);

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
   * @brief Populates a corresponding instance of HaloInfo
   * @param halo the HaloInfo instance to populate
   * @pre halo != NULL
   */
  void GetHaloInfo(HaloInfo *halo);

  /**
   * @brief
   * @param haloParticles
   */
  void GetDIYHaloParticleItemsVector(
          std::vector<DIYHaloParticleItem> &haloParticles);

  /**
   * @brief Registers a DIY data-type to represent static halo data.
   * @param dtype pointer to the DIY data type
   */
 static void CreateDIYHaloType(DIY_Datatype *dtype);

 /**
  * @brief Registers a DIY data-type to represent halo particles.
  * @param dtype pointer to the DIY data type
  */
 static void CreateDIYHaloParticleType(DIY_Datatype *dtype);


  int Count;                    // A count used for book-keeping the number of
                                // times the halo is propagated as a zombie.

  int Tag;                      // The tag/ID of the halo
  int TimeStep;                 // The time-step of this halo
  int HaloType;                 // The type of the halo, see HaloTypeEnum
  int OwnerBlockId;             // The (DIY) global block ID that owns this
                                // halo
  REAL Redshift;                // The corresponding red-shift of the halo
  REAL HaloMass;                // Mass of the halos
  POSVEL_T Center[3];           // The halo-center
  POSVEL_T MeanCenter[3];       // Alternate halo-center definition
  POSVEL_T AverageVelocity[3];  // The average velocity of the halo
  std::set< ID_T > ParticleIds; // The global particle IDs of the halo


private:

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
