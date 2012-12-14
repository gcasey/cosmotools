/**
 * @brief
 */
#ifndef HALOMERGERTREEKERNEL_H_
#define HALOMERGERTREEKERNEL_H_

#include "CosmologyToolsMacros.h"

#include <vector> // For STL vector


namespace cosmotk
{

// Forward declarations
class Halo;

class HaloMergerTreeKernel
{
public:
  HaloMergerTreeKernel();
  virtual ~HaloMergerTreeKernel();

  /**
   * @brief Given two sets of halos at different timesteps, this method
   * will compute the corresponding merger-tree.
   * @param t1 the time-step for the 1st halo set.
   * @param haloSet1 the 1st halo set.
   * @param t2 the time-step for the 2nd halo set.
   * @param haloSet2 the 2nd halo set.
   */
  void ComputeMergerTree(
      const int t1, std::vector< Halo > &haloSet1,
      const int t2, std::vector< Halo > &haloSet2);

  /**
   * @brief Updates the halo-evolution tree, from the current similarity matrix.
   */
  void UpdateHaloEvolutionTree();

protected:
  std::vector< int > HaloSimilarityMatrix; // Flat 2-D matrix
  int MergerTreeThreshold;

private:
  DISABLE_COPY_AND_ASSIGNMENT(HaloMergerTreeKernel);

};

} /* namespace cosmotk */
#endif /* HALOMERGERTREEKERNEL_H_ */
