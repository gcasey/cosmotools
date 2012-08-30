/**
 * @brief A simple test for the CosmologyTools parser.
 */
#include <iostream>
#include <sstream>
#include <mpi.h>

#include "CosmologyToolsMacros.h"
#include "CosmologyToolsConfiguration.h"

/**
 * @brief Program main
 * @param argc the argument counter
 * @param argv the argument vector
 * @return rc test return code rc=0 if succss, else rc!=0.
 */
int main(int argc, char **argv)
{
  int rc = 0;
  int rank = -1;
  int N    = 0;

  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD,&N);

  std::string configFile = std::string(argv[1]);
  std::cout << "Rank=" << rank << " config file: " << configFile << std::endl;
  std::cout.flush();

  cosmotk::CosmologyToolsConfiguration CosmoToolsConfig;
  CosmoToolsConfig.SetCommunicator(MPI_COMM_WORLD);
  CosmoToolsConfig.SetConfigFile(configFile);
  CosmoToolsConfig.ParseFile();

  std::ostringstream oss;
  oss.clear();
  oss << "config-rank-" << rank << ".dat";
  CosmoToolsConfig.WriteConfiguration( oss.str() );

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
  return( rc );
}



