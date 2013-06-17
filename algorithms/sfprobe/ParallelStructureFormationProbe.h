/**
 * @class ParallelStructureFormationProbe
 * @brief A class that extends the functionality implemented by the
 * StructureFormationProbe for MPI distributed particle datasets.
 */
#ifndef PARALLELSTRUCTUREFORMATIONPROBE_H_
#define PARALLELSTRUCTUREFORMATIONPROBE_H_

#include "StructureFormationProbe.h"
#include "CosmoToolsMacros.h"

#include <mpi.h>

namespace cosmologytools {

class ParallelStructureFormationProbe : public StructureFormationProbe
{
public:
  ParallelStructureFormationProbe();
  virtual ~ParallelStructureFormationProbe();

  // In-line functions
  SetMacro(Communicator,MPI_Comm);

  /**
   * @brief Tesselates the grid
   * @see StructureFormationProbe::Tesselate
   */
  virtual void Tesselate();

  /**
   * @brief Barrier synchronization among MPI processes.
   * @post Blocks until all processes hit the barrier.
   */
  void Barrier()
    {MPI_Barrier(this->Communicator);};

  /**
   * @brief Returns the number of ranks associated with the given communicator.
   * @return N the number of ranks
   * @post N >= 1
   */
  int GetNumberOfRanks()
    {
    int size;
    MPI_Comm_size(this->Communicator,&size);
    return( size );
    }

  /**
   * @brief Returns the MPI rank of this process.
   * @return R the rank of this process
   * @post R >= 0 && R < this->GetNumberOfRanks()
   */
  int Rank()
    {
    int rank;
    MPI_Comm_rank(this->Communicator,&rank);
    return rank;
    }

protected:
  MPI_Comm Communicator;

private:
  DISABLE_COPY_AND_ASSIGNMENT(ParallelStructureFormationProbe);
};

} /* namespace hacctools */
#endif /* PARALLELSTRUCTUREFORMATIONPROBE_H_ */
