/**
 * @brief A singleton class that consists of generic functionality used
 * throughout the implementation of the GenericIO readers/writers.
 */
#ifndef GENERICIOUTILITIES_H_
#define GENERICIOUTILITIES_H_

#include "CosmologyToolsMacros.h"
#include "GenericIODefinitions.hpp"

#include <vector>

namespace cosmotk
{

class GenericIOUtilities
{
public:

  /**
   * @brief Checks whether this machine is big endian or not
   * @return status true iff big endian, else false
   */
  static bool IsBigEndian()
    {
    const uint32_t one = 1;
    return !(*((char *)(&one)));
    }

  /**
   * @brief Round-Robins a user-supplied number of blocks to a bucket with the
   * given index, idx. The resulting blocks that are assigned to bucket in
   * query are inserted in the given STL vector.
   * @param bucketIdx the bucket index (in).
   * @param NumBuckets the total number of buckets (in).
   * @param NumBlocks the total number of blocks (in).
   * @param assigned a list of assigned blocks (in/out).
   */
  static void RoundRobin(
      const int bucketIdx, const int NumBuckets,
      const int NumBlocks,
      std::vector<int> &assigned);

  /**
   * @brief Checks if the endian of this machine matches the endian of the
   * file encoded in the global header.
   * @param GH pointer to the global header
   * @return status true if the endian matches, else, false.
   */
  static bool DoesFileEndianMatch(GlobalHeader *GH);

  /**
   * @brief Swaps the endian of the data pointed to by Addr.
   * @param Addr user-supplied address of the data to swap (in/out)
   * @param N the number of bytes of the data to swap (in)
   * @param buf buffer space used to swap data (optional)
   * @note When this method is called within a tight loop it is advised that
   * a buffer is provided and allocated externally. Otherwise, the method will
   * do the necessary buffering internally but will have a cost of allocating
   * and de-allocating small memory for potentially a large number of times.
   */
  static void SwapEndian(void *Addr, const int N, char *buf=NULL);

  /**
   * @brief Swaps the endian of the given global header.
   * @param GH pointer to the global header to swap.
   */
  static void SwapGlobalHeader(GlobalHeader *GH);

  /**
   * @brief Swaps the endian of the Variable header.
   * @param VH pointer to the variable header to swap.
   */
  static void SwapVariableHeader(VariableHeader *VH);

  /**
   * @brief Swaps the endian of the rank header.
   * @param RH pointer to the rank header to swap.
   */
  static void SwapRankHeader(RankHeader *RH);

  /**
   * @brief Given the variable information, it detects the corresponding
   * underlying primitive type to represent the data associated with the
   * variable.
   * @param vinfo the variable information of the variable in query.
   * @return t the primitive type.
   * @see GenericIOPrimitiveTypes.
   */
  static int DetectVariablePrimitiveType(const VariableInfo &vinfo);

  /**
   * @brief Allocates and returns a pointer (void*) to a buffer wherein the
   * data of the corresponding variable can be stored.
   * @param vinfo the variable information of the variable in query.
   * @param numElements the number of elements to allocate.
   * @return bufPtr pointer to allocated buffer.
   * @post bufPtr != NULL.
   */
  static void* AllocateVariableArray(
            const VariableInfo &vinfo,const int numElements);

protected:
  GenericIOUtilities();
  virtual ~GenericIOUtilities();

private:
  DISABLE_COPY_AND_ASSIGNMENT(GenericIOUtilities);
};


} /* namespace cosmotk */
#endif /* GENERICIOUTILITIES_H_ */
