/**
 * @class ForwardHaloTracker
 * @brief A class that implements functionality for in-situ halo tracking. The
 * simulation specifies the frequency, f, at which the halo-finder is invoked.
 * Each time the halo finder is invoked, the resulting mapping of halo tags and
 * global IDS is cached.
 */
#ifndef FORWARDHALOTRACKER_H_
#define FORWARDHALOTRACKER_H_

#include "CosmologyToolsMacros.h"
#include <mpi.h> // For MPI_Comm definition

// C++ includes
#include <set>
#include <vector>
#include <string>

// CosmologyTools includes
#include "Halo.h"
#include "DistributedHaloEvolutionTree.h"
#include "ParallelHaloMergerTree.h"

namespace cosmologytools {

// Forward declarations
class HaloDataInformation;
class TemporalHaloInformation;
class CosmoHaloFinderP;

class ForwardHaloTracker
{
public:
  ForwardHaloTracker();
  virtual ~ForwardHaloTracker();

  /**
   * @brief Registers particles at the given time-step
   * @param tstep the time-step of the particles
   * @param dt the actual time at tstep
   * @param xLoc the x-location of the particles
   * @param yLoc the y-location of the particles
   * @param zLoc the z-location of the particles
   * @param xVel the x-velocity of the particles
   * @param yVel the y-velocity of the particles
   * @param zVel the z-velocity of the particles
   * @param potential the particle potential vector potential for each particle
   * @param id the global IDs of the particles the IDs of each particle
   * @param mask the mask of the particles
   * @param state array used internally by the halofinder
   * @param N the total number of particles
   */
  void RegisterParticles(
      const INTEGER tstep, const REAL redShift,
      POSVEL_T* px, POSVEL_T* py, POSVEL_T *pz,
      POSVEL_T* vx, POSVEL_T* vy, POSVEL_T *vz,
      POSVEL_T* mass, POTENTIAL_T* potential, ID_T* id,
      MASK_T* mask,
      STATUS_T* state,
      INTEGER N);

  // In-line parameters
  GetNSetMacro(BoxLength,REAL);
  GetNSetMacro(NG,INTEGER);
  GetNSetMacro(NDIM,INTEGER);
  GetNSetMacro(PMIN,INTEGER);
  GetNSetMacro(LinkingLength,REAL);
  GetNSetMacro(Communicator,MPI_Comm);
  GetNSetMacro(MergerTreeFileName,std::string);
  GetNSetMacro(MergerTreeThreshold,int);
  GetNSetMacro(MergerTreeFileFormat,int);

  /**
   * @brief Tracks halos.
   * @note This method is intended to be called repeateadly
   */
  void TrackHalos();

  /**
   * @brief Barrier synchronization among all processes on the communicator
   * associated with this instance of the communicator.
   * @note If a communicator is not set, by default, MPI_COMM_WORLD is assumed.
   */
  void Barrier(){ MPI_Barrier(this->Communicator); };

protected:

  /**
   * @brief Extracts the halo information, i.e., the globalIds and
   * corresponding halo tags from an instance of the cosmo halo-finder.
   * @param hinfo the halo information object
   * @param hfinder the halo finder object
   * @pre hinfo != NULL
   * @pre hfinder != NULL
   */
  void GetHaloInformation(
      HaloDataInformation* hinfo, CosmoHaloFinderP* hfinder);

  /**
   * @brief Initializes internal data-structures.
   * @post this->Initialized == true
   */
  void Initialize();

  bool Initialized;

  // Halo finder parameters
  REAL BoxLength;     // length of the box domain
  INTEGER NG;         // size of overlap, percentage of the BoxLength
  INTEGER NDIM;       // number of grid points
  INTEGER PMIN;       // minimum number of halos per particle
  REAL LinkingLength; // linking-length to use for FOF


  INTEGER TimeStep;
  REAL RedShift;

  std::string MergerTreeFileName;
  int MergerTreeThreshold;
  int MergerTreeFileFormat;

  // Registered particles, these pointers are owned by the caller!
  POSVEL_T* Px;
  POSVEL_T* Py;
  POSVEL_T* Pz;
  POSVEL_T* Vx;
  POSVEL_T* Vy;
  POSVEL_T* Vz;
  POSVEL_T* Mass;
  POTENTIAL_T* Potential;
  ID_T* Id;
  MASK_T* Mask;
  STATUS_T* State;
  INTEGER NumberOfParticles;

  MPI_Comm Communicator; // The MPI communicator to use

  TemporalHaloInformation *TemporalHaloData; // Pointer to temporal halo data
  cosmotk::DistributedHaloEvolutionTree *HaloEvolutionTree;
  cosmotk::ParallelHaloMergerTree *HaloMergerTree;

private:

  /**
   * @brief Updates the merger-tree
   */
  void UpdateMergerTree();

  /**
   * @brief Executes the halo-finder.
   * @param haloFinder pointer to a halo-finder instance.
   */
  void ExecuteHaloFinder(CosmoHaloFinderP* haloFinder);

  DISABLE_COPY_AND_ASSIGNMENT(ForwardHaloTracker);
};

} /* namespace cosmogolytools */
#endif /* FORWARDHALOTRACKER_H_ */
