/**
 * @brief ParallelHaloMergerTree extends the HaloMergerTreeKernel and
 * implements functionality for building incrementally a merger-tree
 * on halos that are distributed among several processes.
 */
#ifndef PARALLELHALOMERGERTREE_H_
#define PARALLELHALOMERGERTREE_H_

#include "HaloMergerTreeKernel.h"

#include "CosmologyToolsMacros.h"
#include "DistributedHaloEvolutionTree.h"

// MPI
#include <mpi.h>

// STL data-structures
#include <map>
#include <vector>

namespace cosmotk
{

// Forward declarations
class Halo;
class HaloNeighborExchange;

// Data-structure to store halos based on hash-code
typedef std::map<std::string,Halo> HaloHashMap;

class ParallelHaloMergerTree: public HaloMergerTreeKernel
{
public:
  ParallelHaloMergerTree();
  virtual ~ParallelHaloMergerTree();

  /**
   * @brief Get/Set the communicator
   */
  GetNSetMacro(Communicator,MPI_Comm);

  /**
   * @brief Given two sets of halos at different time-steps, this method
   * updates the merger-tree, represented by a user-supplied halo evolution
   * tree instance.
   * @param t1 the time-step of the first set of halos (in).
   * @param haloSet1 array of halos at t1 (in).
   * @param M the number of halos at t1 (in).
   * @param t2 the time-step of the 2nd set of halos (in).
   * @param haloSet2 array of halos at t2 (in).
   * @param N the number of halos at t2 (in).
   * @param t the halo evolution tree (out).
   * @pre (t1 < t2).
   * @pre (M >= 1) && (N >= 1).
   * @pre (haloSet1 != NULL) && (haloSet2 != NULL).
   * @pre (t != NULL).
   */
  virtual void UpdateMergerTree(
       const int t1, Halo *haloSet1, const int M,
       const int t2, Halo *haloSet2, const int N,
       DistributedHaloEvolutionTree *t);

  /**
   * @brief Computes the global number of birth events detected.
   * @return N the total number of births.
   * @note This method is collective, all ranks must call it.
   */
  int GetTotalNumberOfBirths();

  /**
   * @brief Computes the global number of re-birth events detected.
   * @return N the total number of re-births.
   * @note This method is collective, all ranks must call it.
   */
  int GetTotalNumberOfRebirths();

  /**
   * @brief Computes the global number of merge events detected.
   * @return N the total number of merges.
   * @note This method is collective, all ranks must call it.
   */
  int GetTotalNumberOfMerges();

  /**
   * @brief Computes the global number of split events detected.
   * @return N the total number of splits.
   * @note This method is collective, all ranks must call it.
   */
  int GetTotalNumberOfSplits();

  /**
   * @brief Computes the global number of death events detected.
   * @return N the total number of deaths.
   * @note This method is collective, all ranks must call it.
   */
  int GetTotalNumberOfDeaths();

  /**
   * @brief Barrier synchronization among all processes
   */
  void Barrier() { MPI_Barrier(this->Communicator); };

protected:
   MPI_Comm Communicator;

   // Keeps an incremental count of the number of nodes, i.e., halos per
   // time-step including zombies. Used for calculation of Global IDs.
   ID_T NumberOfNodes;

   // TemporalHalos data-structure is a vector of size 2
   // TemporalHalos[0] holds the list of halos at the previous time-step
   // TemporalHalos[1] holds the list of halos at the current time-step
   // Both lists are global in the sense that they contain, the neighbor
   // halos as well.
   std::vector< std::vector< cosmotk::Halo > > TemporalHalos;

   // Indices used to access the TemporalHalos data-structure
   int CurrentIdx;
   int PreviousIdx;

   HaloNeighborExchange *NeighborExchange;

   /**
    * @brief Assigns a global ID to the given halos.
    * @param halos pointer to the list of halos.
    * @param NumHalos number of halos.
    */
   void AssignGlobalIds(Halo* halos, const int NumHalos);

   /**
    * @brief Updates halo evolution tree accordingly to treat death events.
    * @param t pointer the DistributedHaloEvolutionTree.
    */
   void HandleDeathEvents(DistributedHaloEvolutionTree *t);


private:
   DISABLE_COPY_AND_ASSIGNMENT(ParallelHaloMergerTree);
};

} /* namespace cosmotk */
#endif /* PARALLELHALOMERGERTREE_H_ */
