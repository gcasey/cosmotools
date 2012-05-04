/**
 * @class ExtentPartitioner
 * @brief A class that provides functionality for partitioning extents
 */
#ifndef EXTENTPARTITIONER_H_
#define EXTENTPARTITIONER_H_

#include "CosmologyToolsMacros.h"
#include <vector>

namespace cosmologytools {

class ExtentPartitioner
{
public:
  ExtentPartitioner();
  virtual ~ExtentPartitioner();

  // Inline methods
  SetVector6Macro(GlobalExtent,int);
  GetNSetMacro(NumberOfPartitions,int);

  /**
   * @brief Partitions the global extent to sub-extents
   */
  void Partition();

  /**
   * @brief Returns the extent of the partition corresponding to the given idx.
   * @param idx the index of the partition in query
   * @param ext user-supplied storage array where the extent will be stored.
   * @note the ext stores [imin imax jmin imax kmin kmax]
   * @pre (idx >= 0) && (idx < this->NumberOfPartitions)
   */
  void GetExtent(const int idx, int ext[6]);

protected:
  int GlobalExtent[6];
  int NumberOfPartitions;

  std::vector<int> PartitionExtents;

  /**
   * @brief Returns the total number of partitions
   * @return N the total number of partitions
   */
  int GetTotalNumberOfPartitions()
    {return this->PartitionExtents.size()/6;};

  /**
   * @brief Given an extent, this method computes the longest dimension
   * @param ext the extent in query
   * @return dim the longest dimension
   * @post dim can be either 1,2,3 corresponding to i,j,k respectively.
   */
  int GetLongestDimension(int ext[6]);

  /**
   * @brief Appends the given extent to the end of the PartitionExtents list.
   * @param ext the grid extent added to the list.
   */
  void AppendExtent(int ext[6]);

  /**
   * @brief Replaces the extent at the given location with ext
   * @param idx the index to PartitionExtents list
   * @param ext the ext to be replaced with
   * @pre (idx >= 0) && (idx < this->NumberOfPartitions)
   */
  void ReplaceExtent(const int idx, int ext[6]);

  /**
   * @brief Splits the parent extent in two sub-extents along the given dimension.
   * @param parent the parent extent to be split
   * @param s1 the first sub-extent
   * @param s2 the second sub-extent
   * @param splitDimension the dimension at which the grid will be split.
   */
  void SplitExtent(int parent[6],int s1[6],int s2[6],int splitDimension);

private:
  DISABLE_COPY_AND_ASSIGNMENT(ExtentPartitioner);

};

} /* namespace hacctools */
#endif /* EXTENTPARTITIONER_H_ */
