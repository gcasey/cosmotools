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
void GenericIOMPIReader::Open()
{
  // sanity checks
  assert("pre: empty filename!" && (!this->FileName.empty()) );

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

}

//------------------------------------------------------------------------------
void GenericIOMPIReader::AllocateInternalReaders(const int N)
{
  assert("pre: Internal Readers should be NULL!" &&
		 (this->InternalReaders==NULL) );
  assert("pre: NumberOfReaders(N) > 0" && (N > 0) );

  this->InternalReaders = new GenericIOReader*[N];
  assert("pre: Could not allocate internal readers array!" &&
		  (this->InternalReaders != NULL));

  for( int i=0; i < N; ++i )
    {
	this->InternalReaders[ i ] = new GenericIOMPIReader();
    } // END for all readers
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

      // Read checksum
      uint64_t checksum;
      offSet += bytesize;
      this->Read(&checksum,CRCSize,offSet,"variable checksum");

      // Verify Checksum
      if( !GenericIOUtilities::VerifyChecksum(dataPtr,bytesize,checksum) )
      	{
    	throw std::runtime_error(
    	  "Variable checksum failed for variable " + this->Vars[varIdx].Name);
      	} // END if checksum

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
