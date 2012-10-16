/**
 * @brief A concrete object that handles parsing and storage of the cosmology
 * tools parameters.
 */

#ifndef COSMOLOGYTOOLSCONFIGURATION_H_
#define COSMOLOGYTOOLSCONFIGURATION_H_

#include "CosmologyToolsMacros.h"


// C/C++ includes
#include <map>
#include <string>
#include <vector>

// MPI include
#include <mpi.h>

namespace cosmotk
{

class CosmologyToolsConfiguration
{
public:
  CosmologyToolsConfiguration();
  virtual ~CosmologyToolsConfiguration();

  GetNSetMacro(ConfigFile,std::string);
  GetNSetMacro(Version,double);
  GetNSetMacro(Visualization,bool);
  GetNSetMacro(VisualizationFrequency,int);
  GetNSetMacro(Server,std::string);
  GetNSetMacro(Port,int);
  GetNSetMacro(Communicator,MPI_Comm);

  /**
   * @brief Parses the configuration file
   * @return status true if file is successfully parsed, else false.
   */
  bool ParseFile();

  /**
   * @return N Total number of analysis tools
   */
  int GetNumberOfAnalysisTools();

  /**
   * @brief Get the analysis tool instance name at the given index
   * @param idx the index of the tool
   * @return name the name of the analysis tool
   * @pre idx >= 0 && idx < this->GetNumberOfAnalysisTools()
   * @note This is the internal AnalysisTool name used to determine which
   * class to instantiate.
   */
  std::string GetToolInstanceName(const int idx);

  /**
   * @brief Get the analysis tool name at the given index
   * @param idx the index of the tool
   * @return name the name of the analysis tool
   * @pre idx >= 0 && idx < this->GetNumberOfAnalysisTools()
   */
  std::string GetToolName(const int idx);

  /**
   * @brief Get the parameters of the given tool in a dictionary
   * @param toolName the analysis tool in query
   * @param parameters a key,value pair for the analysis tool parameters
   */
  void GetToolParameters(
        std::string toolName, Dictionary &parameters);

  /**
   * @brief Checks if the tool with the given name exists
   * @param toolName the name of the tool to check
   * @return status true if the tool exists, else, false
   */
  bool ToolExists(std::string toolName);

  /**
   * @brief Checks if the tool with the prescribed name is enabled.
   * @param toolName the name of the tool
   * @return status true if enabled, else false
   */
  bool IsToolEnabled(std::string toolName);

  /**
   * @brief Clears all internal data
   */
  void Clear();

  /**
   * @brief Writes configuration file to disk at the given filename
   * @note Mainly used for debugging
   */
  void WriteConfiguration(std::string file);

protected:
  MPI_Comm Communicator;

  std::string RawConfigData; // The raw configuration data
  std::string ConfigFile; // cosmology tools configuration file
  double Version;        // version of the configuration

  // Visualization parameters
  bool Visualization;         // Indicates whether visualization is enabled
  std::string Server;         // Name or IP address of the server to connect to
  int Port;                   // Port where to connect to visualize outputs
  int VisualizationFrequency; // Frequency at which to update the visualization

  std::vector< std::string >     ToolNames;
  std::map< std::string, bool > ToolStatus;


  std::map< std::string, Dictionary > ToolToDictionary;

  /**
   * @brief Reads in the configuration file to a string
   */
  void ReadInConfigurationFile();

  /**
   * @brief Parses the raw configuration data
   * @pre this->RawConfigData != ""
   */
  bool ParseRawData();

  /**
   * @brief Interprets the given string parameter as a boolean
   * @param s the string parameter to interpret
   * @return true or false
   * @pre s != ""
   */
  bool Parsebool(std::string s);

  /**
   * @brief Registers the analysis tool with the given name and status
   * @param toolName the name of the analysis tool
   * @param status the status (i.e., enabled or disabled)
   */
  void RegisterAnalysisTool(std::string toolName, bool status);

private:
  DISABLE_COPY_AND_ASSIGNMENT(CosmologyToolsConfiguration);
};

} /* namespace cosmotk */
#endif /* COSMOLOGYTOOLSCONFIGURATION_H_ */
