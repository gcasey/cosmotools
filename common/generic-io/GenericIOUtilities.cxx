#include "GenericIOUtilities.h"

#include <cassert>

namespace cosmotk
{

GenericIOUtilities::GenericIOUtilities()
{
  // TODO Auto-generated constructor stub
}

//------------------------------------------------------------------------------
GenericIOUtilities::~GenericIOUtilities()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void GenericIOUtilities::RoundRobin(
    const int bucketIdx, const int NumBuckets,
    const int NumBlocks,
    std::vector<int> &assigned)
{
  assert( "pre: NumBlocks > 0" && (NumBlocks > 0) );
  assert( "pre: NumBuckets > 0" && (NumBuckets > 0) );
  assert( "pre: bucketIdx is out-of bounds" &&
          (bucketIdx >= 0) && (bucketIdx < NumBuckets));

  assigned.clear();

  if(NumBuckets < NumBlocks)
    {

    for(int blkIdx=0; blkIdx < NumBlocks; ++blkIdx)
      {
      if( (blkIdx%NumBuckets) == bucketIdx )
        {
        assigned.push_back(blkIdx);
        } // END if
      } // END for all blocks

    } // END if
  else if( NumBuckets > NumBlocks )
    {

    if( bucketIdx < NumBlocks )
      {
      assigned.push_back(bucketIdx);
      } // END if

    } // END else-if
  else
    {
    assigned.push_back(bucketIdx);
    } // END else
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapEndian(
        void *Addr, const int N, char *buf)
{
  assert("pre: Addr is NULL!" && (Addr != NULL) );
  assert("pre: N > 0" && (N > 0) );

  // STEP 0: Setup buffer. Either the user has supplied a swap buffer
  // externally to use, or we will handle the buffering internally.
  bool bufferedExternally = true;
  if( buf == NULL )
    {
    bufferedExternally = false;
    buf = new char[N];
    }

  // STEP 1: swap bytes
  for(int srcOffSet=N-1, idx=0; srcOffSet >= 0; --srcOffSet, ++idx)
    {
    buf[idx] = *( (char*)Addr+srcOffSet);
    } // END for all bytes

  // STEP 2: Copy data to original
  memcpy(Addr, (void*)buf, N);

  // STEP 3: if buffering internally, clean up
  if( !bufferedExternally )
    {
    delete [] buf;
    }
}

} /* namespace cosmotk */
