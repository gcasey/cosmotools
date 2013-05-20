#include "GenericIOMPIReader.h"

#include "CRC64.h"
#include "GenericIOUtilities.h"
#include "MPIUtilities.h"

// C/C++ includes
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <fstream>
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
  this->IOStrategy   = FileIOMPI;
  this->ProxyEnabled = false;
  this->SplitMode    = false;
  this->InternalReaders = NULL;
}

//------------------------------------------------------------------------------
GenericIOMPIReader::~GenericIOMPIReader()
{
  if( this->SplitMode )
    {
    for(int i=0; i < this->NumberOfFiles; ++i)
      {
      if( this->InternalReaders[i] != NULL)
        {
        delete this->InternalReaders[i];
        }
      }
    delete [] this->InternalReaders;
    } // END if split-mode
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::OpenAndReadHeader( bool skipBlockHeaders )
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

  // STEP 2: Read Global & Variable header
  this->ReadHeader();

  // STEP 3: Index the variables
  this->IndexVariables();

//
//  // STEP 4: Detect file type, i.e., if the data is written to separate files
//  switch( this->Rank )
//    {
//    case 0:
//      {
//      this->DetermineFileType();
//      int split = (this->SplitMode)? 1 : 0;
//      MPI_Bcast(&split,1,MPI_INTEGER,0,this->Communicator);
//      if( this->SplitMode )
//        {
//        skipBlockHeaders = true;
//
//        // Setup to read metadata file in process 0
//        this->AssignedBlocks.push_back( 0 );
//        this->ReadBlockHeaders();
//
//        this->ReadBlockToFileMap();
//
//        // Clear temporary data used for reading the block-to-file mapping
//        this->AssignedBlocks.clear();
//        this->RH.clear();
//        this->VH.clear();
//        } // END if splitmode
//      } // END rank==0
//      break;
//    default:
//      {
//      int split = -1;
//      MPI_Bcast(&split,1,MPI_INTEGER,0,this->Communicator);
//      this->SplitMode = (split==1)? true : false;
//      if( this->SplitMode )
//        {
//        skipBlockHeaders = true;
//        this->ReadBlockToFileMap();
//        } // END if splitmode
//      } // END default
//    } // END switch
//  MPI_Barrier(this->Communicator);
//
//  // STEP 5: Round robin assignment
//  GenericIOUtilities::RoundRobin(
//    this->Rank,this->NumRanks,this->GH.NRanks,this->AssignedBlocks);
//
//  // STEP 6: Setup internal readers, if split mode
//  this->SetupInternalReaders();
//  MPI_Barrier(this->Communicator);
//
//  // STEP 7: Read block headers
//  if( !skipBlockHeaders )
//    {
//    this->ReadBlockHeaders();
//    }
//  MPI_Barrier(this->Communicator);
}

//------------------------------------------------------------------------------
int GenericIOMPIReader::GetNumberOfElements()
{
  int NElements = 0;
  if( this->SplitMode && this->ProxyEnabled )
    {
    for(int i=0; i < this->NumberOfFiles; ++i)
      {
      assert("pre: internal reader is NULL" &&
              (this->InternalReaders[i] != NULL));

      NElements += this->InternalReaders[ i ]->GetNumberOfElements();
      } // END for all files
    } // END if reading in split mode
  else
    {
    unsigned int blkIdx=0;
    for(; blkIdx < this->AssignedBlocks.size(); ++blkIdx)
      {
      NElements += this->GetNumberOfElementsForBlock(blkIdx);
      } // END for all blocks
    }
  return( NElements );
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::SetupInternalReaders()
{
  // STEP 0: Short-circuit if we are reading a single file
  if( !this->SplitMode )
    {
    return;
    }

  this->RH.clear();
  this->VH.clear();

  // STEP 1: Construct all readers, based on number of separate files that
  // we have to read.
  this->InternalReaders = new GenericIOMPIReader *[this->NumberOfFiles];

  std::ostringstream oss;
  for(int i=0; i < this->NumberOfFiles; ++i)
    {
    oss.clear(); oss.str("");
    oss << this->FileName << "#" << i;
    GenericIOMPIReader *iReader = new GenericIOMPIReader();
    iReader->SetCommunicator( this->Communicator );
    iReader->SetFileName(oss.str());
    iReader->OpenAndReadHeader( /* skipBlock headers */ true );
    assert("internal reader should not be in SplitMode" &&
           !iReader->IsSplitMode());

    /* Manually prescribe block assignement for this reader */
    iReader->ClearBlockAssignment();
    this->InternalReaders[ i ] = iReader;
    assert("manually assign blocks to internal reader!" &&
            iReader->GetNumberOfAssignedBlocks()==0);
    }

  // STEP 2: Enable proxying calls to internal readers.
  // Indicate that calls to this class should be proxied to the underlying
  // internal readers.
  this->ProxyEnabled = true;

  // STEP 3: Prescribe block assignment. Note, block assignment is done on
  // the total number of blocks, in OpenAndReadHeader. The total number of
  // blocks is implicitely given by the total number of elements read from
  // the block-to-file map, see ReadBlockToFileMap() method.
  for( unsigned int block=0; block < this->AssignedBlocks.size(); ++block)
    {
    int globalBlockIdx = this->AssignedBlocks[ block ];

    // Sanity checks!
    assert("ERROR: Cannot map block to file!" &&
     (this->BlockToFileMap.find(globalBlockIdx)!=this->BlockToFileMap.end()));
    assert("ERROR: Cannot map block to idx within file!" &&
     (this->BlockToIdxWithinFile.find(globalBlockIdx)!=
      this->BlockToIdxWithinFile.end()));

    int fileIdx = this->BlockToFileMap[ globalBlockIdx ];
    assert("ERROR: fileIdx is out-of-bounds" &&
            (fileIdx >= 0) && (fileIdx < this->NumberOfFiles) );

    // Get block Idx within file
    int idxWithinFile = this->BlockToIdxWithinFile[ globalBlockIdx ];

    this->InternalReaders[ fileIdx ]->AssignBlock( idxWithinFile );
    } // END for all assigned blocks

  // STEP 4: Read the block headers for each reader on each process
  for(int i=0; i < this->NumberOfFiles; ++i)
    {
    this->InternalReaders[ i ]->ReadBlockHeaders();
    }

  // STEP 5: Proxy variable headers to proxy reader
  for(int i=0;i < this->InternalReaders[0]->GetNumberOfVariablesInFile();++i)
    {
    this->VH.push_back(this->InternalReaders[0]->GetVariableHeader(i) );
    } // END for all variables

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
void GenericIOMPIReader::ReadData()
{
  if( this->SplitMode && this->ProxyEnabled )
    {
    this->ReadSplitModeData();
    }
  else
    {
    this->ReadSingleFileData();
    }
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadSplitModeData()
{
  // STEP 0: Propagate variables to internal readers
  int nskip = 0;
  for(int i=0; i < this->NumberOfFiles; ++i)
    {
    assert("pre: Internal reader should not be NULL!" &&
             (this->InternalReaders[ i ] != NULL) );

    this->InternalReaders[ i ]->ClearVariables();

    // Determine offset in to variable array to write
    if( i > 0 )
      {
      nskip += this->InternalReaders[ i-1 ]->GetNumberOfElements();
      }

    for(unsigned int varIdx=0; varIdx < this->Vars.size(); ++varIdx)
      {
      // sanity check!
      assert("pre: Data for variable is NULL!" &&
          (this->Vars[varIdx].Data != NULL) );

      // Get the variable size
      size_t vsize = this->Vars[varIdx].Size;

      // Get pointer where reader "i" will start filling in data for this
      // variable
      void *dataPtr =
          static_cast<char*>(this->Vars[varIdx].Data)+(nskip*vsize);
      assert("pre: dataPtr is NULL!" && (dataPtr != NULL) );

      this->InternalReaders[ i ]->AddVariable(
          this->GetVariableInfo(varIdx),dataPtr);

      } // END for all variables

    // Read all the data from this file
    this->InternalReaders[ i ]->ReadData();
    } // END for all files
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadSingleFileData()
{
  // STEP 0: Get the total number of elements in an array for this rank
  int N = this->GetNumberOfElements();

  // STEP 1: Loop through all registered variables and read them in
  for(unsigned int varIdx=0; varIdx < this->Vars.size(); ++varIdx)
    {
    // sanity check!
    assert( "pre: Data for variable is NULL!" &&
            (this->Vars[varIdx].Data != NULL) );

    // Ensure the variable has extra space
    if(!this->Vars[varIdx].HasExtraSpace)
      {
      throw std::runtime_error(
          "pre: Extra space is required for reading!");
      } // END if extra space

    // Get the variable size
    size_t vsize = this->Vars[varIdx].Size;

    // pointer to data
    void *dataPtr = this->Vars[varIdx].Data;

    // Get the variable index, used to calculate the offset in the file
    int vidx = this->GetVariableIndex( this->Vars[ varIdx ].Name );
    if(vsize != this->VH[vidx].Size)
      {
      std::cerr << "Variable size mismatch for var: "
                << this->Vars[ varIdx ].Name << std::endl;
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

    //
    // TODO: Perform checksum of the variable here
    //

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
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::Close()
{
  MPI_File_close(&this->FH);
  if( this->SplitMode )
    {
    for(int i=0; i < this->NumberOfFiles; ++i )
      {
      this->InternalReaders[ i ]->Close();
      } // END for all files
    } // END if in split mode
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadVariableHeader(
        const int idx, VariableHeader& vh)
{
  uint64_t offSet = this->GH.VarsStart + idx*sizeof(VariableHeader);
  assert("pre: detected variable offset out-of-bounds!" &&
            offSet < this->EntireHeader.size()-CRCSize );

  // Copy the bytes of the variable header from the raw header data
  memcpy(&vh,&this->EntireHeader[offSet],sizeof(VariableHeader));

  if(this->SwapEndian)
    {
    GenericIOUtilities::SwapVariableHeader(&vh);
    }
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadVariableHeaders()
{
  assert( "pre: file has no variables!" && (this->GH.NVars > 0) );

  this->VH.resize( this->GH.NVars );
  for(int i=0; i < this->GH.NVars; ++i )
    {
    this->ReadVariableHeader(i, this->VH[i] );
    } // END for all variables
}

//------------------------------------------------------------------------------
void GenericIOMPIReader::ReadHeader()
{
 // Internal attributes. Each element in the `attribs` array represents a
 // particular attribute, e.g., whether we should swap endian, etc., as
 // indicated below. The reason for storing these attributes in an array
 // instead of individual ints is so that we can send all these attributes
 // with a single broadcast to all ranks.
 int attribs[3]=
   { 0, // indicates whether to swap or not
     0, // indicates the entire header size,including the CRC checksum
     0  // indicates whether an error occured
   };

 // Integers corresponding to indices in the `attribs` array for
 const int SWAP        = 0;
 const int HEADER_SIZE = 1;
 const int ERROR       = 2;

 // Read in attributes
 switch(this->Rank)
   {
   case 0:
     // Read the global header
     this->Read(&this->GH,sizeof(GlobalHeader),0,"GlobalHeader");

     // Read entire header & its checksum
     attribs[HEADER_SIZE] = this->GH.HeaderSize+CRCSize;
     this->EntireHeader.resize(attribs[HEADER_SIZE],0xFE/* poison */);
     this->Read(&this->EntireHeader[0],attribs[HEADER_SIZE],0,"EntireHeader");

     // header checksum -- CRC is endian independent. It must be verified
     // before byte-swapping.
     if(crc64_omp(&this->EntireHeader[0],attribs[HEADER_SIZE])!=(uint64_t)-1)
       {
       attribs[ERROR] = 1;
       }

     // Byte-swap header if necessary
     if( !GenericIOUtilities::DoesFileEndianMatch(&this->GH) )
       {
       this->SwapEndian = true;
       attribs[SWAP]    = 1;
       }
     else
       {
       this->SwapEndian = false;
       attribs[SWAP]    = 0;
       }
     MPI_Bcast(attribs,3,MPI_INTEGER,0,this->Communicator);
     break;
   default:
     MPI_Bcast(attribs,3,MPI_INTEGER,0,this->Communicator);
     this->SwapEndian = (attribs[SWAP]==1)? true : false;
     this->EntireHeader.resize( attribs[HEADER_SIZE] );
   } // END switch

 // Check for errors
 if(attribs[ERROR]==1)
   {
   throw std::runtime_error("Header CRC checksum failed!");
   }

 // Broadcast the raw bytes of the entire header
 assert("pre: headers has not been properly allocated!" &&
           (this->EntireHeader.size()==attribs[HEADER_SIZE]) );
 MPI_Bcast(
   &this->EntireHeader[0],attribs[HEADER_SIZE],MPI_CHAR,0,this->Communicator);

 // Ensure broadcast of the header was successful
 if(crc64_omp(&this->EntireHeader[0],attribs[HEADER_SIZE])!=(uint64_t)-1)
   {
   attribs[ERROR] = 1;
   }
 int errors = 0;
 MPI_Allreduce(
     &attribs[ERROR],&errors,1,MPI_INTEGER,MPI_SUM,this->Communicator);
 if(errors > 0 )
   {
   throw std::runtime_error("Error broadcasting header to ranks!");
   }

 // Extract the global header
 this->GH = *(GlobalHeader*)(&this->EntireHeader[0]);
 if( this->SwapEndian )
   {
   GenericIOUtilities::SwapGlobalHeader(&this->GH);
   }

 // Read the variable headers
 this->ReadVariableHeaders();
 this->Barrier();
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
