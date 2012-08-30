/**
 * @brief CosmologyToolsManager to keep persistent information through the
 * life-time of the simulation.
 */
#ifndef COSMOLOGYTOOLSMANAGER_H_
#define COSMOLOGYTOOLSMANAGER_H_

#include "CosmologyToolsMacros.h"
#include <mpi.h>

#include <cassert>
#include <string>
#include <vector>

namespace cosmotk {

// Forward declaration
class AnalysisTool;
class CosmologyToolsConfiguration;
class SimulationParticles;

class CosmologyToolsManager
{
public:

  /**
   * @brief Default Constructor
   */
  CosmologyToolsManager();

  /**
   * @brief Destructor
   */
  ~CosmologyToolsManager();

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
   */
  void SetDomainParameters(
      REAL boxlength, INTEGER ghostoverlap, INTEGER NDIM);

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
      REAL *mass, POTENTIAL_T *potential,
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

  // Domain parameters
  REAL BoxLength;
  INTEGER GhostOverlap;
  INTEGER NDIM;

  // Data-structure to store the simulation particles
  SimulationParticles *Particles;

  // List of Analysis tools
  std::vector<AnalysisTool*> AnalysisTools;

  /**
   * @brief Clears the analysis tools
   */
  void ClearAnalysisTools();

  /**
   * @brief Parses the configuration file on rank 0 and distributes parameters
   * to all ranks.
   */
  void ParseConfigurationFile();

  /**
   * @brief Creates analysis tools
   */
  void CreateAnalysisTools();

  DISABLE_COPY_AND_ASSIGNMENT(CosmologyToolsManager);
};

} // END cosmotk namespace

#endif
