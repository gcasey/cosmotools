#include "ExtentPartitioner.h"

#include <cassert>
#include <cmath>
#include <queue>

namespace hacctools {

ExtentPartitioner::ExtentPartitioner()
{
  this->NumberOfPartitions  = 2;
}

//-----------------------------------------------------------------------------
ExtentPartitioner::~ExtentPartitioner()
{
  this->PartitionExtents.clear();
}

//-----------------------------------------------------------------------------
void ExtentPartitioner::Partition()
{
  this->PartitionExtents.clear();

  this->AppendExtent(this->GlobalExtent);
  std::queue<int> wrkQueue;
  wrkQueue.push( 0 );

  int ex[6];
  int s1[6];
  int s2[6];

  while( this->GetTotalNumberOfPartitions() < this->NumberOfPartitions )
    {
    int idx = wrkQueue.front();
    wrkQueue.pop();

    this->GetExtent( idx, ex );
    int ldim = this->GetLongestDimension( ex );
    this->SplitExtent( ex, s1, s2, ldim );
    this->ReplaceExtent( idx, s1);
    this->AppendExtent(s2);
    wrkQueue.push( idx );
    wrkQueue.push(this->GetTotalNumberOfPartitions()-1);
    }

}

//-----------------------------------------------------------------------------
int ExtentPartitioner::GetLongestDimension(int ext[6])
{
  int ilength = (ext[1]-ext[0])+1;
  int jlength = (ext[3]-ext[2])+1;
  int klength = (ext[5]-ext[4])+1;

  if ((ilength >= jlength) && (ilength >= klength))
    {
    return 1;
    }
  else if ((jlength >= ilength) && (jlength >= klength))
    {
    return 2;
    }
  else if ((klength >= ilength) && (klength >= jlength))
    {
    return 3;
    }

  assert( "pre: could not find longest dimension" && false );
  return 0;
}

//-----------------------------------------------------------------------------
void ExtentPartitioner::GetExtent(const int idx, int ext[6])
{
  assert("pre: idx is out-of-bounds!" &&
      (idx >= 0) && (idx <this->NumberOfPartitions) );

  for( int i=0; i < 6; ++i )
    {
    ext[ i ] = this->PartitionExtents[idx*6+i];
    }
}

//-----------------------------------------------------------------------------
void ExtentPartitioner::AppendExtent(int ext[6] )
{
  for( int i=0; i < 6; ++i )
    {
    this->PartitionExtents.push_back( ext[ i ] );
    }
}

//-----------------------------------------------------------------------------
void ExtentPartitioner::ReplaceExtent(const int idx, int ext[6])
{
  assert("pre: idx is out-of-bounds!" &&
        (idx >= 0) && (idx < this->NumberOfPartitions) );

  for( int i=0; i < 6; ++i )
    {
    this->PartitionExtents[idx*6+i]=ext[i];
    }
}

//-----------------------------------------------------------------------------
void ExtentPartitioner::SplitExtent(
    int parent[6], int s1[6], int s2[6],int splitDimension)
{
  int numNodes = 0;
  int mid      = -1;
  int minIdx   = -1;
  int maxIdx   = -1;

  for( int i=0; i < 6; ++i )
    {
    s1[ i ] = s2[ i ] = parent[ i ];
    }

  switch( splitDimension )
    {
    case 1:
      minIdx = 0;
      maxIdx = 1;
      break;
    case 2:
      minIdx = 2;
      maxIdx = 3;
      break;
    case 3:
      minIdx = 4;
      maxIdx = 5;
      break;
    default:
      assert( "Cannot split extent: Undefined split dimension!" && false);
    }

  numNodes      = (parent[maxIdx]-parent[minIdx]) + 1;
  mid           = static_cast<int>(std::floor(0.5*numNodes));
  s1[ maxIdx ]  = (mid < s1[minIdx])? (s1[minIdx]+mid) : mid;
  s2[ minIdx ]  = (mid < s1[minIdx])? (s1[minIdx]+mid) : mid;
}

} /* namespace hacctools */
