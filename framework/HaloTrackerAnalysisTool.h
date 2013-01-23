/**
 * @brief A concrete instance of AnalysisTool that implements functionality for
 * running the halo tracker within the in situ framework. The halo-tracker
 * manages an incremental parallel merger-tree internally that is being
 * updated as the simulation progresses.
 */
#ifndef HALOTRACKERANALYSISTOOL_H_
#define HALOTRACKERANALYSISTOOL_H_

#include "CosmologyToolsMacros.h"
#include "AnalysisTool.h"
#include "ForwardHaloTracker.h"

namespace cosmotk
{

class HaloTrackerAnalysisTool : public AnalysisTool
{
public:

  /**
   * @brief Default constructor
   */
  HaloTrackerAnalysisTool();

  /**
   * @brief Destructor
   */
  virtual ~HaloTrackerAnalysisTool();

  /**
   * @brief Parses the parameters of this analysis tool
   */
  virtual void ParseParameters();

  /**
   * @brief Executes the halo-tracker on the given particle input
   * @param particles pointer to the simulation particles
   * @pre particles != NULL
   * @pre this->Communicator != MPI_COMM_NULL
   */
  virtual void Execute(SimulationParticles *particles);

  /**
   * @brief Writes the output of the halo-finder
   * @see AnalysisTool::WriteOutput
   */
  virtual void WriteOutput();

  /**
   * @brief Returns the information of this AnalysisTool instance
   * @return s a string consisting of information for this class instance.
   */
  virtual std::string GetInformation();

protected:

  // HaloFinder parameters
  REAL LinkingLength;
  INTEGER PMIN;
  INTEGER CenterFinderMethod;
  bool ComputSODHalos;
  REAL RHO_C;
  REAL INITIAL_SOD_MASS;
  REAL MIN_RADIUS_FACTOR;
  REAL MAX_RADIUS_FACTOR;
  INTEGER NUMBER_OF_BINS;
  INTEGER FOF_SIZE_THRESHOLD;

  INTEGER MERGER_TREE_FILE_FORMAT;
  INTEGER MERGER_TREE_THRESHOLD;
  std::string MERGER_TREE_FILE;

  // The underlying halo-tracker algorithm used
  cosmologytools::ForwardHaloTracker *HaloTracker;


private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloTrackerAnalysisTool);
};

} /* namespace cosmotk */
#endif /* HALOTRACKERANALYSISTOOL_H_ */
