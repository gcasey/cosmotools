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
