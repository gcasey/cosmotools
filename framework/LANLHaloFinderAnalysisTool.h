/**
 * @brief A concrete instance of AnalysisTool that implements functionality for
 * running the LANL halofinder within the in situ framework.
 */
#ifndef LANLHALOFINDERANALYSISTOOL_H_
#define LANLHALOFINDERANALYSISTOOL_H_

#include "CosmologyToolsMacros.h"
#include "AnalysisTool.h"
#include "CosmoHaloFinderP.h"

namespace cosmotk
{

enum CenterFinder
{
  AVERAGE = 0,
  CENTER_OF_MASS = 1,
  MBP = 2,
  MCP = 3,

  NUMBER_OF_CENTER_FINDER_METHODS
};

class LANLHaloFinderAnalysisTool : public AnalysisTool
{
public:

  /**
   * @brief Default constructor
   */
  LANLHaloFinderAnalysisTool();

  /**
   * @brief Destructor
   */
  virtual ~LANLHaloFinderAnalysisTool();

  /**
   * @brief Parses the parameters of this analysis tool
   */
  virtual void ParseParameters();

  /**
   * @brief Executes the halo-finder on the given particle input
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
   * @brief Returns the information of this
   * @return
   */
  virtual std::string GetInformation();

protected:

  // LANLHaloFinder parameters from the configuration file
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

  cosmologytools::CosmoHaloFinderP *HaloFinder;
private:
  DISABLE_COPY_AND_ASSIGNMENT(LANLHaloFinderAnalysisTool);
};

} /* namespace cosmotk */
#endif /* LANLHALOFINDERANALYSISTOOL_H_ */
