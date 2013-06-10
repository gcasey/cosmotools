/**
 * @brief CosmologyToolsManager to keep persistent information through the
 * life-time of the simulation.
 */
#ifndef InSituAnalysisManager_H_
#define InSituAnalysisManager_H_

#include "CosmoToolsMacros.h"
#include <mpi.h>

#include <cassert>
#include <map>
#include <string>
#include <vector>

namespace cosmotk {

// Forward declaration
class InSituAlgorithm;
class InSituAnalysisConfig;
class SimulationParticles;

class InSituAnalysisManager
{
public:

  /**
   * @brief Default Constructor
   */
  InSituAnalysisManager();

  /**
   * @brief Destructor
   */
  ~InSituAnalysisManager();

  /**
   * @brief Initializes the cosmology environment
   * @param comm user-supplied communicator, default uses MPI_COMM_WORLD
   * @post this->Communicator != MPI_COMM_NULL
   */
  void Initialize(MPI_Comm comm=MPI_COMM_WORLD);

  /**
   * @brief Sets the configuration file to be used by this instance of the
   * cosmology tools.
   * @param configfile configuration file
   */
  void SetAnalysisConfiguration(char *configfile);

  /**
   * @brief Sets the domain parameters for this instance of cosmology tools.
   * @param boxlength the length of the box domain
   * @param ghostoverlap the overlap of ghost particles between partitions and
   * across periodic boundaries
   * @param NDIM the number of points along each dimension in lagrangian space
   * @param XYZPeriodic flag that indicates if the domain is XYZ periodic
   */
  void SetDomainParameters(
      REAL boxlength, INTEGER ghostoverlap, INTEGER NDIM,
      bool XYZPeriodic);

  /**
   * @brief Check if the given time-step is an execution time-step.
   * @param tstep the current time-step
   * @param redshift the redshift at the given time-step
   * @return status true iff one or more analysis tools must execute in this
   * time-step, otherwise, false.
   */
  bool IsExecutionTimeStep(INTEGER tstep, REAL redshift);

  /**
   * @brief Sets the particle information at the given time-step.
   * @param tstep the current time-step
   * @param redshift the redshift at the given time-step
   * @param px the x-coordinate of the particle position vector
   * @param py the y-coordinate of the particle position vector
   * @param pz the z-coordinate of the particle position vector
   * @param vx the x-coordinate of the particle velocity vector
   * @param vy the y-coordinate of the particle velocity vector
   * @param vz the z-coordinate of the particle velocity vector
   * @param mass array of particle masses
   * @param potential array of particle potential
   * @param tags array of particle tags
   * @param mask array of particle masking
   * @param status array of particle status
   * @param N the total number of particles for the given process
   */
  void SetParticles(
      INTEGER tstep, REAL redShift,
      POSVEL_T *px, POSVEL_T *py, POSVEL_T *pz,
      POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz,
      POSVEL_T *mass, POTENTIAL_T *potential,
      ID_T *tags, MASK_T *mask, STATUS_T *status,
      ID_T N);

  /**
   * @brief Executes analysis tasks according to the configuration file.
   * @pre this->Communicator != MPI_COMM_NULL
   */
  void CoProcess();

  /**
   * @brief Performs barrier synchronization.
   */
  void Barrier()
    {
    assert("pre: communicator is NULL!"&&(this->Communicator!=MPI_COMM_NULL));
    MPI_Barrier(this->Communicator);
    }

  /**
   * @brief Finalizes cosmology tools environment
   */
  void Finalize();

private:
  std::string ConfigurationFile;
  MPI_Comm Communicator;
  int Rank;
  int NumRanks;

  // Flag that indicates whether the configuration has been read at the
  // supplied time-step
  bool ConfigurationIsRead;

  // Domain parameters
  REAL BoxLength;
  INTEGER GhostOverlap;
  INTEGER NDIM;
  bool XYZPeriodic;

  // Storage for configuration parameters
  InSituAnalysisConfig *Configuration;

  // Data-structure to store the simulation particles
  SimulationParticles *Particles;

  // List of Analysis tools by name
  std::map<std::string,InSituAlgorithm*> InSituAlgorithms;

  // Timers & other stats
  std::map<std::string,double>      Timers;       // Discrete times at each tstep
  std::map<std::string, INTEGER > EventCounter; // Number of event occurences
  std::map<std::string,double>     GlobalTimers; // Total times (accumulative)

  /**
   * @brief Checks if the communicator associated with this instance of the
   * InSituAnalysisManager is a cartesian communicator.
   * @return status true if the communicator is cartesian, else, false.
   */
  bool IsCartesianCommunicator();

  /**
   * @brief Computes the linear rank of the process with the given position
   * @param i the ith position of the rank in query
   * @param j the jth position of the rank in query
   * @param k the kth position of the rank in query
   * @return r the linear rank.
   * @pre The communicator topology must be cartesian
   */
  int GetRankByPosition(int i, int j, int k);

  /**
   * @brief Given the cartesian position of this rank, this method computes
   * the neighbors of this rank, including periodic neighbors.
   * @param pos the position of this rank (in)
   * @param neighbor the ranks of the 26 neighbors.
   * @note The domain is XYZ-periodic, hence, each rank will have exactly 26
   * neighbors, i.e., 6 face neighbors, 12 edge neighbors, 8 corner neighbors.
   * @pre Assumes a cartesian communicator.
   */
  void ComputeRankNeighbors(
          int pos[3], int neighbor[26]);

  /**
   * Computes this rank's block origin, i.e., min, and size based on the
   * decomposition size and the position.
   * @param decompSize the cartesian decomposition dimensions (in)
   * @param pos the position of this rank within the cartesian communicator (in)
   * @param min the origin (x,y,z) coordinates of this block (out)
   * @param size the size of this block (out)
   */
  void GetBlockBounds(
       int decompSize[3], int pos[3], float min[3], float size[3]);

  /**
   * @brief Describe the decomposition for DIY.
   */
  void SetupDIYDecomposition();

  /**
   * @brief Gather send vector on all processes
   * @param send the send of vector on this process
   * @param rcv the rcv vector consisting of the data from all processes
   * @note A helper method to collect statistics from all processes
   */
  void GatherVector(
      std::vector<double> &send, std::vector<double> &rcv);

  /**
   * @brief Starts a timer for the given event at the given timestep
   * @param tstep the current timestep
   * @param eventName the event name
   */
  void StartTimer(INTEGER tstep, std::string eventName);

  /**
   * @brief Ends the timer for the givent event at the given timestep
   * @param tstep the current timestep
   * @param eventName the event name
   */
  void EndTimer(INTEGER tstep, std::string eventName);

  /**
   * @brief Gathers all timers from all processes to rank 0 and write them
   * in a file.
   */
  void FinalizeTimers();

  /**
   * @brief Clears the analysis tools
   */
  void ClearInSituAlgorithms();

  /**
   * @brief Parses the configuration file on rank 0 and distributes parameters
   * to all ranks.
   */
  void ParseConfigurationFile();

  /**
   * @brief Creates analysis tools
   */
  void CreateInSituAlgorithms();

  /**
   * @brief Checks if the tool exists in the list of analysis for this
   * InSituAnalysisManager instance.
   * @param toolName the name of the tool in query
   * @return status true if the
   */
  bool ToolExists(const std::string &toolName);

  /**
   * @brief Get the analysis tool instance for the tool in query.
   * @param toolName the name of the tool in query.
   * @return tool pointer to analysis tool instance
   * @post tool != NULL
   */
  InSituAlgorithm* GetToolByName(const std::string &toolName);

  DISABLE_COPY_AND_ASSIGNMENT(InSituAnalysisManager);
};

} // END cosmotk namespace

#endif
