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
   * @param potential the particle potential ?
   * @param id the global IDs of the particles ?
   * @param mask the mask of the particles ?
   * @param state ?
   * @param N the total number of particles
   */
  void RegisterParticles(
      const int tstep, const double redShift,
      REAL* px, REAL* py, REAL *pz,
      REAL* vx, REAL* vy, REAL *vz,
      REAL* mass, REAL* potential, INTEGER* id,
      INTEGER* mask, INTEGER* state,
      INTEGER N);

  /**
   * @brief Explicitly set the tracker time-steps when the tracker is invoked.
   * @param tsteps array of time-steps
   * @param N the number of time-steps
   */
  void SetExplicitTrackerTimeSteps(INTEGER *tsteps, const INTEGER N);

  /**
   * @brief Checks if the supplied timesteps is a tracker timestep. If explicit
   * timesteps are provided, then the method checks if the given timestep is
   * in the list of pre-scribed time-steps that the tracker should execute.
   * Otherwise, the frequency is used to determine if the tracker should run.
   * @param tstep the current time-step
   * @return status true if it is a tracker time-step, else false.
   */
  bool IsTrackerTimeStep(const int tstep);


  // In-line halo-finder parameters
  GetNSetMacro(PMin,INTEGER);
  GetNSetMacro(LinkingLength,REAL);
  GetNSetMacro(RL,REAL);
  GetNSetMacro(Overlap,int);

  // In-line parameter for the tracker
  GetNSetMacro(Frequency,int);
  GetNSetMacro(Communicator,MPI_Comm);

  /**
   * @brief Updates the merger tree according to the user-supplied frequency
   * @param tstep the discrete time-step
   * @note if tstep%this->Frequency != 0 the method returns immediately.
   */
  void UpdateMergerTree(const int tstep);

  /**
   * @brief Barrier synchronization among all processes on the communicator
   * associated with this instance of the communicator.
   * @note If a communicator is not set, by default, MPI_COMM_WORLD is assumed.
   */
  void Barrier(){ MPI_Barrier(this->Communicator); };

protected:

  /**
   * @brief Extracts the halo information, i.e., the globalIds and corresponding
   * halo tags from an instance of the cosmo halo-finder.
   * @param hinfo the halo information object
   * @param hfinder the halo finder object
   * @pre hinfo != NULL
   * @pre hfinder != NULL
   */
  void GetHaloInformation(
      HaloDataInformation* hinfo, CosmoHaloFinderP* hfinder);

  // Halo finder parameters
  INTEGER PMin;       // minimum number of halos per particle
  REAL LinkingLength; // linking-length to use for FOF
  REAL RL;            // The physical box size (i.e., domain size)
  INTEGER Overlap;    // the ghost overlap

  // Registered particles, these pointers are owned by the caller!
  std::vector< REAL > Px;
  std::vector< REAL > Py;
  std::vector< REAL > Pz;
  std::vector< REAL > Vx;
  std::vector< REAL > Vy;
  std::vector< REAL > Vz;
  std::vector< REAL > Mass;
  std::vector< REAL > Potential;
  std::vector< INTEGER > Id;
  std::vector< short unsigned int > Mask;
  std::vector< INTEGER > State;
  INTEGER NumberOfParticles;

  // Tracker parameters
  bool UseExplicitTimeSteps;
  std::set< INTEGER > TimeSteps;

  MPI_Comm Communicator; // The MPI communicator to use
  int Frequency;        // The execution frequency of the halo-finder/tracker.

  TemporalHaloInformation *TemporalHaloData;
private:
  DISABLE_COPY_AND_ASSIGNMENT(ForwardHaloTracker);
};

} /* namespace cosmogolytools */
#endif /* FORWARDHALOTRACKER_H_ */
