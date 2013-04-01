/**
 * @brief Merger-tree kernel that supports incremental building of a
 * halo-evolution tree instance, used to capture the halo formation
 * history during the course of a simulation.
 */
#ifndef HALOMERGERTREEKERNEL_H_
#define HALOMERGERTREEKERNEL_H_

// CosmologyTools Macros
#include "CosmologyToolsMacros.h"
#include "DistributedHaloEvolutionTree.h" // For DistributedHaloEvolutionTree
#include "Halo.h" // For Halo

// C/C++ includes
#include <cassert> // For assert()
#include <vector>  // For STL vector
#include <set>     // For STL set


namespace cosmotk
{

class HaloMergerTreeKernel
{
public:
  HaloMergerTreeKernel();
  virtual ~HaloMergerTreeKernel();

  /**
   * @brief Get/Set to turn on/off verbosity. Default is true.
   */
  GetNSetMacro(Verbose,bool);

  /**
   * @brief Get/Set macro for the merger-tree. Default is 10.
   */
  GetNSetMacro(MergerTreeThreshold,int);

  /**
   * @brief Get/Set macro for the zombie cut-off. Default is 5.
   */
  GetNSetMacro(ZombieCutOff,int);

  /**
   * @return Get number of births that occurred for the current set of halos
   * spanning two timesteps.
   */
  GetMacro(NumberOfBirths,int);

  /**
   * @return Get number of re-births that occurred for the current set of halos
   * spanning two timesteps.
   */
  GetMacro(NumberOfRebirths,int);

  /**
   * @brief Returns the number of dead halos found
   * @return N the number of dead halos
   * @note Dead halos are halos that existed in a previous
   * timestep but, siezed to exist in the current timestep.
   */
  int GetNumberOfDeadHalos()
    {return this->DeadHalos.size(); };

  /**
   * @brief Returns the number of split halos found
   * @return N the number of halos that were split
   * @note Split halos are halos that in the curr
   */
  int GetNumberOfSplitHalos()
    {return this->SplitHalos.size(); };

  /**
   * @brief Returns the number of merger-halos.
   * @return N the number of merger halos.
   * @note Mergers are halos that have resulted from merging two or
   * more halos from the previous timestep.
   */
  int GetNumberOfMergeHalos()
    {return this->MergeHalos.size(); };

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

protected:

  int Timesteps[2]; // The two time-steps of halo-information.
  int Sizes[2];     // The number of halos for each of the two time-steps.
  Halo *Halos1;     // The set of halos at time-step t1.
  Halo *Halos2;     // The set of halos at time-step t2.

  std::set< int > DeadHalos; // index to the Halos1 array of halos that do
                             // not exist in Halos2

  std::set< int > SplitHalos; // index to the Halos1 array for each halo that
                              // was split in two (or more) halos in the
                              // current timestep.

  std::set< int > MergeHalos; // index to the Halos2 array for each halo at the
                              // current timestep that resulted from a merge
                              // of two (or more) halos in the previous
                              // timestep.

  int NumberOfBirths;         // Counts the number of births that occured.

  int NumberOfRebirths;       // Counts the number of re-births.

  // Container for columns in the similarity matrix that need no further
  // processing. For example, when a birth is detected, this column no longer
  // needs to be checked.
  std::set< int > ProcessedColumns;

  // The MatrixColumnSum stores the sum of each column in the similarity
  // matrix. Based on this sum, we can infer, if a new halo is "born" or
  // if we have a merge event.
  std::vector< int > MatrixColumnSum;

  // The MatrixRowSum stores the sum of each row. Based on this sum, we
  // can infer, if a death occured, if we have a split.
  std::vector< int > MatrixRowSum;

  // The similarity matrix stores the percent similarity of halo at a previous
  // timestep, i.e., Halos1, and a halo at the current timestep, i.e., the
  // array Halos2.
  std::vector< int > HaloSimilarityMatrix; // Flat 2-D matrix

  bool Verbose;     // Parameter that indicate whether to print out debug
                     // output and other information.

  int ZombieCutOff; // Parameter that indicates after how many times of
                     // propagating a halo as a zombie, the halo should
                     // be designated as a zombie and bypass any checks.

  // The user-supplied percent threshold, default is set to 50%
  int MergerTreeThreshold;

  /**
   * @brief Inserts the halo into the merger-tree
   * @param haloPtr pointer to the halo being inserted
   * @param bitmask mask encoding event type
   * @param t pointer to the merger-tree
   * @pre haloPtr != NULL
   * @pre t != NULL
   * @pre !t->HasNode( haloPtr->GetHashCode() )
   * @post t->HashNode( haloPtr->GetHashCode() )
   */
  void InsertHalo(
      Halo *haloPtr, unsigned char bitmask,
      DistributedHaloEvolutionTree *t)
    {
    assert("pre: NULL halo pointer!" && (haloPtr != NULL) );
    assert("pre: NULL merger-tree pointer!" && (t != NULL) );
    assert("pre: node already exists!" &&
            (!t->HasNode(haloPtr->GetHashCode())));

    HaloInfo hinfo;
    haloPtr->GetHaloInfo( &hinfo);
    t->InsertNode( hinfo, bitmask );

    assert("post: inserted node, not found in the tree!" &&
            (t->HasNode(haloPtr->GetHashCode())) );
    };

  /**
   * @brief Checks if the given ID corresponds to an index in the similarity
   * matrix, corresponding to a halo that has been designated as dead.
   * @param idx the index of the halo in query.
   * @return status true if the the corresponding halo has been designated as
   * dead, else, false.
   */
  bool IsHaloDead(const int idx)
    {return(this->DeadHalos.find(idx) != this->DeadHalos.end());};

  /**
   * @brief Checks if the halo corresponding to the given index is a merger
   * @param idx the index of the halo in query.
   * @return status true if the corresponding halo has been designated as
   * merger, else, false.
   */
  bool IsHaloMerge(const int idx)
    {return(this->MergeHalos.find(idx) != this->MergeHalos.end());};

  /**
   * @brief Checks if the halo correponding to the given index has been split
   * @param idx the index of the halo in query.
   * @return status true if the correponding halo is split, else false.
   */
  bool IsHaloSplit(const int idx)
    {return(this->SplitHalos.find(idx) != this->SplitHalos.end());};

  /**
   * @brief Applies majority rule to determine if two halos should be linked.
   * @param percentOverlap the percent of overlap w.r.t. the previous halo.
   * @return status true if the halos match else false.
   */
  bool MajorityRuleCheck(const int percentOverlap)
    {return(percentOverlap >= this->MergerTreeThreshold);};

  /**
   * @brief Registers the halos at the given timesteps.
   * @param t1 the time-step of the first halo set.
   * @param haloSet1 array of halos at t1.
   * @param M the number of halos in haloSet1.
   * @param t2 the time-step of the second halo set.
   * @param haloSet2 array of halos at t2.
   * @param N the number of halos in haloSet2.
   * @pre (t1 < t2).
   * @pre (M >= 1) && (N >= 1).
   * @pre (haloSet1 != NULL) && (haloSet2 != NULL).
   */
  void RegisterHalos(
       const int t1, Halo *haloSet1, const int M,
       const int t2, Halo *haloSet2, const int N);

  /**
   * @brief Given two sets of halos at different time-steps, this method
   * will compute the corresponding merger-tree matrix. The merger-tree is
   * defined in an MxN matrix, H, wherein H(i,j) > MergerThreeThreshold
   * indicates that halo_i \in t_1 is a parent of halo_j \in t_2.
   */
  void ComputeMergerTree();

  /**
   * @brief Processes the similarity matrix and updates the merger-tree.
   * @param t merger-tree being updated.
   * @see DetectEvent
   */
  void UpdateHaloEvolutionTree(DistributedHaloEvolutionTree *t);

  /**
   * @brief Detects the event for a given halo in a previous timestep and
   * updates the merger-tree accordingly.
   * @param row the row-index in the similarity matrix corresponding to the
   * halo from a previous timestep which is being processed.
   * @param rowEvent the event that was detected based on
   * @param prevHalo
   * @param t
   */
  void DetectEvent(
      const int row, const int rowEvent,
      Halo *prevHalo, DistributedHaloEvolutionTree *t);

  /**
   * @brief Prints the current instance of the similarity matrix.
   * @note Used primarily for debugging.
   */
  void PrintMatrix();

private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloMergerTreeKernel);

};

} /* namespace cosmotk */
#endif /* HALOMERGERTREEKERNEL_H_ */
