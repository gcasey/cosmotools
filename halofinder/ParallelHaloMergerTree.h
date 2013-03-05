/**
 * @brief ParallelHaloMergerTree extends the HaloMergerTreeKernel and
 * implements functionality for building incrementally a merger-tree
 * on halos that are distributed among several processes.
 */
#ifndef PARALLELHALOMERGERTREE_H_
#define PARALLELHALOMERGERTREE_H_

#include "HaloMergerTreeKernel.h"

#include "CosmologyToolsMacros.h"
#include <mpi.h>

// STL data-structures
#include <map>
#include <vector>

namespace cosmotk
{

// Forward declarations
class Halo;

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
    * @brief Barrier synchronization among all processes
    */
   void Barrier() { MPI_Barrier(this->Communicator); };

protected:
   MPI_Comm Communicator;

   // TemporalHalos data-structure is a vector of size 2
   // TemporalHalos[0] holds the list of halos at the previous time-step
   // TemporalHalos[1] holds the list of halos at the current time-step
   // Both lists are global in the sense that they contain, the neighbor
   // halos as well.
   std::vector< std::vector< cosmotk::Halo > > TemporalHalos;

   // Indices used to access the TemporalHalos data-structure
   int CurrentIdx;
   int PreviousIdx;

   /**
    * @brief Prints the similarity matrix on this process.
    * @note Used mostly for debugging.
    */
   void PrintMatrix();

   /**
    * @brief Exchanges halos of this process with neighbor processes
    * @param halos the list of halos in this process
    * @param N the number of halos in this process
    * @param gloabalHalos global halo list consisting of both local and
    * remote halos.
    */
   void ExchangeHalos(
       cosmotk::Halo* halos, const int N,
       std::vector<cosmotk::Halo>& gloabalHalos);

   /**
    * @brief Exchange halo static information, e.g., center, velocity, etc.
    * @param halos list of local halos in this process to send
    * @param N the number of local halos
    * @param haloHash hash to store neighbor halos by halo hash code.
    * @see ExchangeHalos
    */
   void ExchangeHaloInfo(
      cosmotk::Halo* halos, const int N, HaloHashMap &haloHash);

   /**
    * @brief Exchange halo particle IDs for each halo
    * @param halos list of local halos
    * @param N the number of local halos
    * @param haloHash hash of neighbor halos
    * @pre haloHash should not be empty, for each neighboring halo, the hash
    * must consist of the static halo information
    * @see ExchangeHaloInfo, ExchangeHalos
    */
   void ExchangeHaloParticles(
       cosmotk::Halo* halos, const int N, HaloHashMap &haloHash);


private:
   DISABLE_COPY_AND_ASSIGNMENT(ParallelHaloMergerTree);
};

} /* namespace cosmotk */
#endif /* PARALLELHALOMERGERTREE_H_ */
