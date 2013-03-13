#include "GenericIOMPIReader.h"
#include "GenericIOUtilities.h"

// C/C++ includes
#include <cassert>
#include <cstddef>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>

namespace cosmotk
{

//------------------------------------------------------------------------------
GenericIOMPIReader::GenericIOMPIReader()
{
  this->Communicator = MPI_COMM_NULL;
  this->SwapEndian   = false;
  this->Rank         = 0;
  this->NumRanks     = 0;
  this->VH.resize( 0 );
  this->RH.resize( 0 );
  this->AssignedBlocks.resize( 0 );
  this->IOStrategy  = FileIOMPI;
}

//------------------------------------------------------------------------------
GenericIOMPIReader::~GenericIOMPIReader()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::OpenAndReadHeader()
{
  // sanity checks
  assert( "pre: No communicator is supplied" &&
          (this->Communicator != MPI_COMM_NULL) );
  assert("pre: FileName is empty!" && (!this->FileName.empty()) );


  // STEP 0: Get Rank and NumRanks information
  MPI_Comm_rank(this->Communicator, &this->Rank);
  MPI_Comm_size(this->Communicator, &this->NumRanks);

  // STEP 1: Open file, if not successful, throw runtime_error exception
  int rc = MPI_File_open(
                this->Communicator,
                const_cast<char*>(this->FileName.c_str()),
                MPI_MODE_RDONLY,
                MPI_INFO_NULL,
                &this->FH);
  if( rc != MPI_SUCCESS )
    {
    throw std::runtime_error( "Unable to open file: " + this->FileName );
    }

  // STEP 2: Read Global &Variable header
  this->ReadHeader();

  // STEP 3: Index the variables
  this->IndexVariables();

  // STEP 4: Detect file type, i.e., if the data is written to separate files
  switch( this->Rank )
    {
    case 0:
      {
      this->DetermineFileType();
      int split = (this->SplitMode)? 1 : 0;
      MPI_Bcast(&split,1,MPI_INTEGER,0,this->Communicator);
      if( this->SplitMode )
        {
        this->ReadBlockToFileMap();
        } // END if splitmode
      } // END rank==0
      break;
    default:
      {
      int split = -1;
      MPI_Bcast(&split,1,MPI_INTEGER,0,this->Communicator);
      this->SplitMode = (split==1)? true : false;
      if( this->SplitMode )
        {
        this->ReadBlockToFileMap();
        } // END if splitmode
      } // END default
    } // END switch


  // STEP 4: Round robin assignment
  GenericIOUtilities::RoundRobin(
      this->Rank,this->NumRanks,this->GH.NRanks,this->AssignedBlocks);

  // STEP 5: Read block headers
  this->ReadBlockHeaders();
  MPI_Barrier(this->Communicator);
}

//------------------------------------------------------------------------------
int GenericIOMPIReader::GetNumberOfElements()
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
void GenericIOMPIReader::ReadBlockToFileMap()
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
      break;
    default:
      MPI_Bcast(&NumElements,1,MPI_INTEGER,0,this->Communicator);
      rank = new int[NumElements];
      part = new int[NumElements];
      MPI_Bcast(rank,NumElements,MPI_INTEGER,0,this->Communicator);
      MPI_Bcast(part,NumElements,MPI_INTEGER,0,this->Communicator);
    } // END switch

  // STEP 1: All processes
  std::set<int > parts;
  for(int i=0; i < NumElements; ++i)
    {
    parts.insert(part[i]);
    this->BlockToFileMap[ rank[i] ]=part[i];
    }
  this->NumberOfFiles = parts.size();

  // STEP 2: Delete dynamically allocated data
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
void GenericIOMPIReader::ReadData()
{
  // STEP 0: Get the total number of elements in an array for this rank
  int N = this->GetNumberOfElements();

  // STEP 1: Loop through all registered variables and read them in
  for(unsigned int varIdx=0; varIdx < this->Vars.size(); ++varIdx)
    {
    // sanity check!
    assert( "pre: Data for variable is NULL!" &&
            (this->Vars[varIdx].Data != NULL) );

    // Get the variable size
    size_t vsize = this->Vars[varIdx].Size;

    // pointer to data
    void *dataPtr = this->Vars[varIdx].Data;

    // Get the variable index, used to calculate the offset in the file
    int vidx = this->GetVariableIndex( this->Vars[ varIdx ].Name );
    if(vsize != this->VH[vidx].Size)
      {
      throw std::runtime_error(
          "Variable size mismatch for " + this->Vars[ varIdx ].Name);
      }
    assert("pre: cannot find variable index!" &&
            (vidx >= 0) && (vidx < this->VH.size()) );

    // Loop through all blocks and read corresponding variable
    for(unsigned int block=0; block < this->AssignedBlocks.size(); ++block)
      {
      // Get the number of elements in the block
      int NBlockElements = this->GetNumberOfElementsForBlock(block);

      // Calculate the offset in the file for the given variable
      uint64_t offSet = this->GetVariableOffSet(vidx,block);

      // Compute number of bytes to read
      size_t bytesize = NBlockElements*vsize;

      // Read in the block data
      this->Read(dataPtr,bytesize,offSet,this->Vars[varIdx].Name);
      dataPtr = static_cast<char*>(dataPtr)+bytesize;
      } // END for all assigned blocks

    // Swap endian if necessary
    if(this->SwapEndian)
      {
      std::vector<char> swapBuffer;
      swapBuffer.resize(vsize);
      void *ptr = this->Vars[varIdx].Data;
      for(int i=0; i < N; ++i, ptr=static_cast<char*>(ptr)+vsize)
        {
        GenericIOUtilities::SwapEndian(ptr,vsize,&swapBuffer[0]);
        }
      } // END if swap endian

    } // END for all variables

  MPI_Barrier(this->Communicator);
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::Close()
{
  MPI_File_close(&this->FH);
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadVariableHeader(
        const int idx, VariableHeader *vh)
{
  assert("pre: vh!=NULL" && (vh != NULL) );

  std::ostringstream oss;
  oss << "Reading variable " << idx;
  uint64_t offSet = this->GH.VarsStart + idx*sizeof(VariableHeader);
  this->Read(vh,sizeof(VariableHeader),offSet, oss.str() );
  if(this->SwapEndian)
    {
    GenericIOUtilities::SwapVariableHeader(vh);
    }
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadVariableHeaders()
{
  assert( "pre: file has no variables!" && (this->GH.NVars > 0) );

  this->VH.resize( this->GH.NVars );
  for(int i=0; i < this->GH.NVars; ++i )
    {
    this->ReadVariableHeader(i, &this->VH[i] );
    } // END for all variables
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadHeader()
{
 int shouldSwapInt; // int used to broadcast should swap.
 int N = 0;         // the number of variables in the file.
 switch(this->Rank)
   {
   case 0:
     // Read header
     this->Read(&this->GH,sizeof(GlobalHeader),0,"GlobalHeader");

     // Byte-swap header if necessary
     if( !GenericIOUtilities::DoesFileEndianMatch(&this->GH) )
       {
       this->SwapEndian = true;
       shouldSwapInt    = 1;
       GenericIOUtilities::SwapGlobalHeader(&this->GH);
       }
     else
       {
       this->SwapEndian = false;
       shouldSwapInt    = 0;
       }
     this->ReadVariableHeaders();
     N = this->VH.size();
     MPI_Bcast(&shouldSwapInt,1,MPI_INTEGER,0,this->Communicator);
     MPI_Bcast(&N,1,MPI_INTEGER,0,this->Communicator);
     break;
   default:
     MPI_Bcast(&shouldSwapInt,1,MPI_INTEGER,0,this->Communicator);
     this->SwapEndian = (shouldSwapInt==1)? true : false;
     MPI_Bcast(&N,1,MPI_INTEGER,0,this->Communicator);
     this->VH.resize( N );
   } // END switch

   // Broadcast header
   MPI_Bcast(&this->GH,sizeof(GlobalHeader),MPI_BYTE,0,this->Communicator);

   // Broadcast variable header
   MPI_Bcast(
      &this->VH[0],N*sizeof(VariableHeader),MPI_BYTE,0,this->Communicator);
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadBlockHeader(
        const int blkIdx, RankHeader *blockHeader)
{
  assert("pre: blockHeader != NULL" && (blockHeader != NULL) );

  std::ostringstream oss;
  oss << "Reading block header for block " << blkIdx;
  uint64_t offSet = this->GH.RanksStart + blkIdx*sizeof(RankHeader);
  this->Read(blockHeader,sizeof(RankHeader),offSet,oss.str());
  if(this->SwapEndian)
    {
    GenericIOUtilities::SwapRankHeader(blockHeader);
    }
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadBlockHeaders()
{
  this->RH.resize( this->AssignedBlocks.size() );
  for(unsigned int blk=0; blk < this->AssignedBlocks.size(); ++blk)
    {
    int blkIdx = this->AssignedBlocks[ blk ];
    this->ReadBlockHeader(blkIdx,&this->RH[blk]);
    } // END for all assigned blocks
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::Read(
        void *buf, size_t count, off_t offset, const std::string &varName)
{
  while( count > 0 )
    {
    MPI_Status status;
    int rc = MPI_File_read_at(this->FH,offset,buf,count,MPI_BYTE,&status);
    if( rc != MPI_SUCCESS )
      {
      throw std::runtime_error(
          "Unable to read " + varName + " form file " + this->FileName);
      }

    int scount;
    MPI_Get_count(&status,MPI_BYTE,&scount);
    count -= scount;
    buf = ((char*) buf) + scount;
    offset += scount;
    } // END while
}

} /* namespace cosmotk */
