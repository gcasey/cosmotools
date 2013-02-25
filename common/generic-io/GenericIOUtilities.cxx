#include "GenericIOUtilities.h"

#include <iostream>
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
void GenericIOUtilities::SwapGlobalHeader(GlobalHeader *GH)
{
  assert("pre: GH != NULL" && (GH != NULL) );

  std::vector<char> swapBuffer;
  swapBuffer.resize( sizeof(char) );
  for(int i=0; i < MagicSize; ++i )
    {
    GenericIOUtilities::SwapEndian(
        &GH->Magic[i],sizeof(char),&swapBuffer[0]);
    }

  swapBuffer.resize( sizeof(uint64_t) );

  GenericIOUtilities::SwapEndian(
      &GH->HeaderSize,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->NElems,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->Dims[0],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->Dims[1],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->Dims[2],sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->NVars,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->VarsSize,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->VarsStart,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->NRanks,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->RanksSize,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->RanksStart,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->GlobalHeaderSize,sizeof(uint64_t),&swapBuffer[0]);

  swapBuffer.resize( sizeof(double) );

  GenericIOUtilities::SwapEndian(
      &GH->PhysOrigin[0],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysOrigin[1],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysOrigin[2],sizeof(double),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &GH->PhysScale[0],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysScale[0],sizeof(double),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &GH->PhysScale[0],sizeof(double),&swapBuffer[0]);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapVariableHeader(VariableHeader *VH)
{
  assert("pre: VH != NULL" && (VH != NULL) );
  std::vector<char> swapBuffer;
  swapBuffer.resize(sizeof(char));

  for(int i=0; i < NameSize; ++i)
    {
    GenericIOUtilities::SwapEndian(
        &VH->Name[i],sizeof(char),&swapBuffer[0]);
    }

  swapBuffer.resize(sizeof(uint64_t));
  GenericIOUtilities::SwapEndian(
      &VH->Flags,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
        &VH->Size,sizeof(uint64_t),&swapBuffer[0]);
}

//------------------------------------------------------------------------------
void GenericIOUtilities::SwapRankHeader(RankHeader *RH)
{
  assert("pre: RH != NULL" && (RH != NULL) );

  std::vector<char> swapBuffer;
  swapBuffer.resize( sizeof(uint64_t) );

  GenericIOUtilities::SwapEndian(
      &RH->Coords[0],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &RH->Coords[1],sizeof(uint64_t),&swapBuffer[0]);
  GenericIOUtilities::SwapEndian(
      &RH->Coords[2],sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &RH->NElems,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &RH->Start,sizeof(uint64_t),&swapBuffer[0]);

  GenericIOUtilities::SwapEndian(
      &RH->GlobalRank,sizeof(uint64_t),&swapBuffer[0]);
}


//------------------------------------------------------------------------------
bool GenericIOUtilities::DoesFileEndianMatch(GlobalHeader *GH)
{
  assert("pre: GH != NULL" && (GH != NULL)  );

  const char *Magic =
       GenericIOUtilities::IsBigEndian() ? MagicBE : MagicLE;
  const char *MagicInv =
       GenericIOUtilities::IsBigEndian() ? MagicLE : MagicBE;

  std::string magicString =
     std::string(GH->Magic,GH->Magic+MagicSize-1);

  if( magicString != Magic)
   {
   if( magicString == MagicInv )
     {
     return false;
     } // END if swap
   else
     {
     std::cerr << "ERROR: Could not detect file endian!\n";
     } // END else
   } // END if endian does not match file endian
  return true;
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
