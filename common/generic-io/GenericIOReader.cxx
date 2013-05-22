#include "GenericIOReader.h"

#include <iostream>
#include <set>

namespace cosmotk
{

//------------------------------------------------------------------------------
GenericIOReader::GenericIOReader()
{
  this->SwapEndian    = false;
  this->SplitMode     = false;
  this->NumberOfFiles = 1;
}

//------------------------------------------------------------------------------
GenericIOReader::~GenericIOReader()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void GenericIOReader::DetermineFileType()
{
  if(  this->HasVariable("$rank")          &&
       this->HasVariable("$partition"))
    {
    this->SplitMode = true;
    }
}

//------------------------------------------------------------------------------
void GenericIOReader::IndexVariables()
{
  std::string name;
  for(unsigned int i=0; i < this->VH.size(); ++i)
    {
    name = std::string( this->VH[i].Name );
    this->VariableName2IndexMap[ name ]= i;
    }
}

//------------------------------------------------------------------------------
int GenericIOReader::GetVariableIndex(const std::string name)
{
  int idx = -1;
  if( this->VariableName2IndexMap.find(name) !=
       this->VariableName2IndexMap.end())
    {
    idx = this->VariableName2IndexMap[name];
    }
  return( idx );
}

//------------------------------------------------------------------------------
void GenericIOReader::AssignBlock(const int blkIdx)
{
  assert("pre: block index is out of bounds" &&
          (blkIdx >= 0) && (blkIdx < this->GH.NRanks) );
  this->AssignedBlocks.push_back( blkIdx );
}

//------------------------------------------------------------------------------
VariableInfo GenericIOReader::GetFileVariableInfo(const int i)
{
  assert("pre: variable index is out-of-bounds!" &&
          (i >= 0) && (i < this->VH.size() ) );

  VariableInfo VI(
      std::string( this->VH[ i ].Name ),
      this->VH[ i ].Size,
      static_cast<bool>(this->VH[i].Flags & FloatValue),
      static_cast<bool>(this->VH[i].Flags & SignedValue),
      static_cast<bool>(this->VH[i].Flags & ValueIsPhysCoordX),
      static_cast<bool>(this->VH[i].Flags & ValueIsPhysCoordY),
      static_cast<bool>(this->VH[i].Flags & ValueIsPhysCoordZ),
      static_cast<bool>(this->VH[i].Flags & ValueMaybePhysGhost)
      );
  return( VI );
}

//------------------------------------------------------------------------------
void GenericIOReader::ReadBlockToFileMap()
{
  assert("pre: file must be in SplitMode" && this->SplitMode);

  // STEP 0: Rank 0 reads block-to-file mapping and distributes it to all ranks
  int NumElements=0;
  int *rank = NULL;
  int *part = NULL;
  switch(this->Rank)
    {
    case 0:
      NumElements = this->GetNumberOfElements();
      MPI_Bcast(&NumElements,1,MPI_INTEGER,0,this->Communicator);
      rank = new int[NumElements];
      part = new int[NumElements];
      this->AddVariable("$rank",rank);
      this->AddVariable("$partition",part);
      this->ReadData();
      MPI_Bcast(rank,NumElements,MPI_INTEGER,0,this->Communicator);
      MPI_Bcast(part,NumElements,MPI_INTEGER,0,this->Communicator);
      this->ClearVariables();

      // Correct the global header
      this->GH.NRanks = NumElements;
      break;
    default:
      MPI_Bcast(&NumElements,1,MPI_INTEGER,0,this->Communicator);
      rank = new int[NumElements];
      part = new int[NumElements];
      MPI_Bcast(rank,NumElements,MPI_INTEGER,0,this->Communicator);
      MPI_Bcast(part,NumElements,MPI_INTEGER,0,this->Communicator);

      // Correct the global header
      this->GH.NRanks = NumElements;
    } // END switch

  // STEP 1: All processes
  std::set<int> parts;
  std::map<int,int> counter;
  int blkIdxInFile = -1;
  for(int i=0; i < NumElements; ++i)
    {
    parts.insert(part[i]);
    this->BlockToFileMap[ rank[i] ] = part[i];

    blkIdxInFile = -1;
    if( counter.find(part[i]) != counter.end() )
      {
      blkIdxInFile    = counter[ part[i] ];
      counter[ part[i] ] += 1;
      } // END if
    else
      {
      blkIdxInFile = 0;
      counter[ part[i] ] = 1;
      } // END else
    this->BlockToIdxWithinFile[ rank[i] ] = blkIdxInFile;
    } // END for all elements
  this->NumberOfFiles = parts.size();

  // STEP 2: Delete dynamically allocated data
  counter.clear();
  parts.clear();
  if( rank != NULL )
    {
    delete [] rank;
    }
  if( part != NULL )
    {
    delete [] part;
    }
}

//------------------------------------------------------------------------------
int GenericIOReader::GetNumberOfElements()
{
  int N = 0;
  unsigned int blkIdx=0;
  for(; blkIdx < this->AssignedBlocks.size(); ++blkIdx)
    {
    N += this->GetNumberOfElementsForBlock(blkIdx);
    } // END for all blocks
  return( N );
}

//------------------------------------------------------------------------------
int GenericIOReader::GetNumberOfElementsForBlock(const int blkidx)
{
  assert("pre: blkIdx is out-of-bounds!" &&
          (blkidx >= 0) && (blkidx < this->AssignedBlocks.size() ) );
  assert("pre: blkIdx is out-of-bounds!" &&
          (blkidx >= 0) && (blkidx < this->RH.size() ) );
  return( this->RH[blkidx].NElems);
}

//------------------------------------------------------------------------------
uint64_t GenericIOReader::GetVariableOffSet(int vidx, int localBlkIdx)
{
  // Sanity checks!
  assert("pre: variable index is out-of-bounds!" &&
       (vidx >= 0) && (vidx < this->GH.NVars) );
  assert("pre: local block index is out-of-bounds!" &&
       (localBlkIdx >= 0) && (localBlkIdx < this->RH.size() ) );


  int N = this->GetNumberOfElementsForBlock(localBlkIdx);
  uint64_t offSet = this->RH[ localBlkIdx ].Start;
  for(int i=0; i < vidx; ++i)
    {
    offSet += N*this->VH[i].Size + CRCSize;
    }
  return( offSet );
}

} /* namespace cosmotk */
