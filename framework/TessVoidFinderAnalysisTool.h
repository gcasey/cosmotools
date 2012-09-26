/**
 * @brief A concrete instance of AnalysisTool that implements functionality for
 * running the tess void finder within the in situ framework.
 */

#ifndef TESSVOIDFINDERANALYSISTOOL_H_
#define TESSVOIDFINDERANALYSISTOOL_H_

#include "CosmologyToolsMacros.h"
#include "AnalysisTool.h"

#include <mpi.h>

namespace cosmotk
{

class TessVoidFinderAnalysisTool : public AnalysisTool
{
public:
  TessVoidFinderAnalysisTool();
  virtual ~TessVoidFinderAnalysisTool();

  /**
   * @brief Parses the parameters of this analysis tool
   * @see AnalysisTool::ParseParameters
   */
  virtual void ParseParameters();

  /**
   * @brief Executes tess on the given particle input
   * @param particles pointer to the simulation particles
   * @pre particles != NULL
   * @pre this->Communicator != MPI_COMM_NULL
   */
  virtual void Execute(SimulationParticles *particles);

  /**
   * @brief Writes output of this analysis tool
   */
  virtual void WriteOutput();

  /**
   * @brief Returns the information of this AnalysisTool instance
   * @return s information of this tool instance in a string
   */
  virtual std::string GetInformation();

protected:
  REAL GhostFactor;
  REAL MinVol;
  REAL MaxVol;
  REAL CellSize;

  bool Initialized; // Indicates whether tess is initialized
  double *TimeStatistics; // Used by tess profiling

  /**
   * @brief Computes the linear rank of the process with the given position
   * @param i the ith position of the rank in query
   * @param j the jth position of the rank in query
   * @param k the kth position of the rank in query
   * @return r the linear rank.
   * @pre The communicator topology must be cartesian
   */
  int GetRankByPosition(int i, int j, int k)
    {
    int ijk[3];
    ijk[0]=i; ijk[1]=j; ijk[2]=k;
    int rank;
    MPI_Cart_rank(this->Communicator,ijk,&rank);
    return( rank );
    }

  /**
   * @brief Given the cartesian position of this rank, this method computes the
   * neighbors of this rank, including periodic neighbors
   * @param pos the cartesian position of this rank (in)
   * @param neighbors the neighbors of this rank (out)
   * @note Because the domain is XYZ periodic each rank will have exactly
   * 26 neighbors,i.e.,6 face neighbors,12 edge neighbors,8 corner neighbors
   * @pre The communicator topology must be cartesian
   */
  void ComputeRankNeighbors(int pos[3],int neighbors[26]);

  /**
   * @brief Given the cartesian topology of this rank, this method computes
   * the local block bounds.
   * @param decompSize the dimensions of the cartesian topology (in)
   * @param pos the position of this rank within the cartesian communicator (in)
   * @param min the minimum cartesian coordinates of this block (out)
   * @param size the size of this block (out)
   */
  void GetBlockBounds(
      int decompSize[3], int pos[3], float min[3], float size[3]);

  /**
   * @brief Initializes Tess
   * @post this->Initialized == true
   */
  void InitializeTess(SimulationParticles *particles);

  /**
   * @brief Packages the particles positions in a flat array for tess
   * @param particles the simulation particles at the given timestep (in)
   * @param positions the flat array (out)
   * @post positions != NULL
   * @note The particle positions is allocated via malloc, b/c tess is a C
   * library and reallocs internally for its particle exchange. To clear the
   * allocated positions, ClearParticlePositions must be called.
   */
  void PackageParticlePositions(
      SimulationParticles *particles, float ***positions);

  /**
   * @brief Clears the particle positions
   * @param positions the positions vector
   * @pre the positions vector must be constructed via PackageParticlePositions
   */
  void ClearParticlePositions( float **positions );

private:
  DISABLE_COPY_AND_ASSIGNMENT(TessVoidFinderAnalysisTool);
};

} /* namespace cosmotk */
#endif /* TESSVOIDFINDERANALYSISTOOL_H_ */
