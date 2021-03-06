/**
 * @brief An abstract base class that defines the basic properties of an
 * InSituAlgorithm instance.
 */
#ifndef InSituAlgorithm_H_
#define InSituAlgorithm_H_

#include "CosmoToolsMacros.h"

// C/C++ includes
#include <cassert>
#include <set>
#include <string>

// MPI include
#include <mpi.h>

#include "SimulationParticles.h"

namespace cosmotk
{


enum InSituAlgorithmFrequency
{
  EXPLICIT = 0,
  IMPLICIT = 1
};

class InSituAlgorithm
{
public:

  /**
   * @brief Default constructor
   */
  InSituAlgorithm();

  /**
   * @brief Destructor
   */
  virtual ~InSituAlgorithm();

  /**
   * @brief Set/Get the name of this analysis tool instance
   */
  GetNSetMacro(Name,std::string);

  /**
   * @brief Set whether output is generated by this InSituAlgorithm instance
   */
  GetNSetMacro(GenerateOutput,bool);

  /**
   * @brief Set the output file name where the output of this InSituAlgorithm is
   * going to be written.
   */
  GetNSetMacro(OutputFile,std::string);

  /**
   * @brief Set the frequency type, i.e., implicit or explicit.
   */
  GetNSetMacro(FrequencyType,int);

  /**
   * @brief Set the implicit frequency at which this algorithm will be executed.
   */
  GetNSetMacro(ImplicitFrequency,int);

  /**
   * @brief Get the parameters of this InSituAlgorithm instance.
   */
  GetMacro(Parameters,Dictionary);

  /**
   * @brief Set MPI Communicator
   */
  GetNSetMacro(Communicator,MPI_Comm);

  /**
   * @brief Enable/disable this tool instance
   */
  GetNSetMacro(Enabled,bool);
  inline bool IsEnabled() {return this->Enabled;};
  std::string GetStringStatus()
   {return( (this->IsEnabled())? "ENABLED" : "DISABLED");}

  /**
   * @brief Set visibility status, i.e., whether, the algorithm will be visible.
   */
  GetNSetMacro(VisibilityStatus,bool);
  inline bool IsVisible() {return this->VisibilityStatus;};

  /**
   * @brief Sets the parameters to this InSituAlgorithm instance.
   * @param parameters dictionary with user-supplied parameters
   */
  void SetParameters(Dictionary parameters)
    {
    this->Parameters = parameters;
    this->ParseParameters();
    }

  /**
   * @brief Set the domain parameters
   * @param boxLength the length of the box
   * @param NG the size of the ghost-overlap zone
   * @param NDIM number of points in each dimension
   * @param Periodic indicates whether the domain is XYZ periodic or not
   * @note domain is assumed to be square
   */
  void SetDomainParameters(
      REAL boxLength, INTEGER NG, INTEGER NDIM, bool periodic=true)
      {
      this->BoxLength=boxLength;
      this->NG=NG;
      this->NDIM=NDIM;
      this->Periodic=periodic;
      };

  /**
   * @brief Sets the explicit timesteps at which this algorithm will be
   * executed.
   * @param timeSteps array of timesteps
   * @param N the number of timesteps.
   */
  void SetExplicitTimeSteps(INTEGER *timeSteps, int N);

  /**
   * @brief Determines if the algorithm should execute at the current timestep
   * @param ts the current timestep
   * @return status true if the algorithm must execute, else, false.
   */
  bool ShouldExecute(INTEGER ts);

  /**
   * @return rank of process w.r.t the supplied communicator
   */
  int Rank()
    {
    assert("pre: communicator is NULL!"&&(this->Communicator!=MPI_COMM_NULL));
    int rank = 0;
    MPI_Comm_rank(this->Communicator,&rank);
    return( rank );
    }

  /**
   * @brief Performs barrier synchronization.
   */
  void Barrier()
    {
    assert("pre: communicator is NULL!"&&(this->Communicator!=MPI_COMM_NULL));
    MPI_Barrier(this->Communicator);
    }

  /**
   * @brief Parses the parameters from the given dictionary.
   */
  virtual void ParseParameters() = 0;

  /**
   * @brief Executes the analysis tool. Implementation is defined by
   * concrete classes.
   */
  virtual void Execute(SimulationParticles *particles) = 0;

  /**
   * @brief Writes the output to a file. Implementation is defined by
   * concrete classes.
   */
  virtual void WriteOutput() = 0;

  /**
   * @brief Gets the analysis tool information in a string
   * @note Used mainly for debugging.
   * @return the information string
   */
  virtual std::string GetInformation() = 0;

protected:

  // Domain parameters
  REAL BoxLength;
  INTEGER NG;
  INTEGER NDIM;
  bool Periodic;

  // Name of the analysis tool, set in the constructor of concrete classes
  std::string Name;

  // Common analysis tool parameters, read from the configuration file
  bool GenerateOutput;
  std::string OutputFile;

  int FrequencyType;
  int ImplicitFrequency;
  std::set<INTEGER> ExplicitTimeSteps;

  bool Enabled;

  bool VisibilityStatus;

  // Storage of parameters
  Dictionary Parameters;

  // MPI communicator used
  MPI_Comm Communicator;

  /**
   * @brief Returns the value of the parameter with the given key
   * @param key the name of the parameter in query
   * @return the value of the parameter
   */
  double GetDoubleParameter(std::string key);

  /**
   * @brief Returns the value of the parameter with the given key
   * @param key the name of the parameter in query
   * @return the value of the parameter
   */
  int GetIntParameter(std::string key);

  /**
   * @brief Returns the value of the parameter with the given key
   * @param key the name of the parameter in query
   * @return the value of the parameter
   */
  std::set<int> GetIntListParameter(std::string key);

  /**
   * @brief Returns the value of the parameter with the given key
   * @param key the name of the parameter in query
   * @return the value of the parameter
   */
  bool GetBooleanParameter(std::string key);

  /**
   * @brief Returns the value of the parameter with the given key
   * @param key the name of the parameter in query
   * @return the value of the parameter
   */
  std::string GetStringParameter(std::string key);

  /**
   * @brief Parse basic/common parameters for all analysis tools
   * @note This method is intended as a helper method for concrete instances
   * @pre Parameters.empty() == false
   */
  void ParseBasicParameters();

  /**
   * @brief Returns a string representation of the common parameters among all
   * analysis tools.
   * @return str a string representation of the analysis tool basic info
   */
  std::string GetBasicInformation();

private:
  DISABLE_COPY_AND_ASSIGNMENT(InSituAlgorithm);
};

} /* namespace cosmotk */
#endif /* ANALYSISTOOL_H_ */
