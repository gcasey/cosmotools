/**
 * @brief Merger-tree kernel that supports incremental building of a
 * halo-evolution tree instance, used to capture the halo formation
 * history during the course of a simulation.
 */
#ifndef HALOMERGERTREEKERNEL_H_
#define HALOMERGERTREEKERNEL_H_

#include "CosmologyToolsMacros.h"

#include <vector> // For STL vector


namespace cosmotk
{

// Forward declarations
class Halo;
class DistributedHaloEvolutionTree;

class HaloMergerTreeKernel
{
public:
  HaloMergerTreeKernel();
  virtual ~HaloMergerTreeKernel();

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

  std::vector< int > HaloSimilarityMatrix; // Flat 2-D matrix
  int MergerTreeThreshold;


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
   * will compute the corresponding merger-tree. The merger-tree is defined
   * in an MxN matrix, H, wherein H(i,j) > 0 indicates that halo_i \in t_1 is
   * a parent of halo_j \in t_2.
   */
  void ComputeMergerTree();

  /**
   * @brief Updates the halo-evolution tree, from the current similarity matrix.
   * @param t the halo evolution tree instance to update.
   */
  void UpdateHaloEvolutionTree(
          DistributedHaloEvolutionTree *t);

private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloMergerTreeKernel);

};

} /* namespace cosmotk */
#endif /* HALOMERGERTREEKERNEL_H_ */
