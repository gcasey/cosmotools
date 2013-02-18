/**
 * @brief A singleton class that consists of generic functionality used
 * throughout the implementation of the GenericIO readers/writers.
 */
#ifndef GENERICIOUTILITIES_H_
#define GENERICIOUTILITIES_H_

#include "CosmologyToolsMacros.h"

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

protected:
  GenericIOUtilities();
  virtual ~GenericIOUtilities();

private:
  DISABLE_COPY_AND_ASSIGNMENT(GenericIOUtilities);
};


} /* namespace cosmotk */
#endif /* GENERICIOUTILITIES_H_ */
