/**
 * @brief A simple program that opens a file written in the GenericIO format
 * and prints out the number of variables, their corresponding sizes and
 * names.
 */

// C++ includes
#include <cassert>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// MPI include
#include <mpi.h>

// CosmologyTools includes
#include "GenericIODefinitions.hpp"
#include "GenericIOMPIReader.h"
#include "GenericIOPosixReader.h"
#include "GenericIOReader.h"
#include "GenericIOUtilities.h"
#include "MPIUtilities.h"

//==============================================================================
// Global variables
//==============================================================================
int rank;
int size;
MPI_Comm comm = MPI_COMM_WORLD;

/* Command line arguments */
bool OnlyMetaData = true;
bool UsePosix	  = false;
std::string File  = "";

struct DataInformation
{
  cosmotk::VariableInfo VariableInformation;
  void *Data;
};

//==============================================================================
// Global Methods
//==============================================================================
template<class T>
std::string PrintData(void *dataPtr, const int idx)
{
  T* castedPtr = static_cast< T* >( dataPtr );
  std::ostringstream sbuf;
  sbuf.clear();
  sbuf.str("");

  sbuf << std::scientific
            << std::setprecision( std::numeric_limits<T>::digits10 )
            << castedPtr[ idx ] << "\t";
  return( sbuf.str() );
}

//==============================================================================
void ParseArguments(int argc, char **argv)
{
  for(int i=1; i < argc; ++i)
    {
    if(strcmp(argv[i],"--file")==0)
      {
      File = std::string(argv[++i]);
      }
    else if(strcmp(argv[i],"--use-posix")==0)
      {
    UsePosix = true;
      }
    else if(strcmp(argv[i],"--all")==0)
      {
    OnlyMetaData = false;
      }
    } // END for all arguments

  assert("pre: specify [--file <file>]" && (File != "") );
}

/**
 * @brief Program main
 * @param argc argument counter
 * @param argv argument vector
 * @return rc return code
 */
int main(int argc, char **argv)
{
  // STEP 0: Initialize
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(comm,&rank);
  MPI_Comm_size(comm,&size);

  // STEP 1: Get file name to read
  ParseArguments(argc,argv);

  // STEP 2: Read and print vars
  cosmotk::GenericIOReader *reader = NULL;
  if( UsePosix )
    {
  reader = new cosmotk::GenericIOPosixReader();
    }
  else
    {
  reader = new cosmotk::GenericIOMPIReader();
    }

  reader->SetCommunicator(comm);
  reader->SetFileName(File);
  reader->OpenAndReadHeader();

  std::vector< DataInformation > DataVector;
  if( !OnlyMetaData )
  {
  DataVector.resize( reader->GetNumberOfVariablesInFile() );
  }

  int N = reader->GetNumberOfElements();
  cosmotk::MPIUtilities::Printf(
      comm,"# NUMBER OF VARIABLES=%d\n",reader->GetNumberOfVariablesInFile());
  cosmotk::MPIUtilities::SynchronizedPrintf(
      comm,"# NUMBER OF ELEMENTS=%d\n",N);
  MPI_Barrier(comm);

  std::ostringstream oss;
  oss.str("");
  oss.clear();
  for(int i=0; i < reader->GetNumberOfVariablesInFile(); ++i)
    {
  if( !OnlyMetaData )
    {
      DataVector[ i ].VariableInformation = reader->GetFileVariableInfo( i );
    }

    int type = cosmotk::GenericIOUtilities::DetectVariablePrimitiveType(
                              reader->GetFileVariableInfo(i));

    assert( "pre: undefined data type!" &&
            (type >= 0) && (type < cosmotk::NUM_PRIMITIVE_TYPES) );

    oss.str("");
    oss << reader->GetVariableName(i) << " ("
        << cosmotk::PRIMITIVE_NAME[type] << ");";
    cosmotk::MPIUtilities::Printf(comm,"# VARIABLE=%s\n",oss.str().c_str());

    if( !OnlyMetaData )
    {
    DataVector[ i ].Data =
    cosmotk::GenericIOUtilities::AllocateVariableArray(
        DataVector[ i ].VariableInformation,N);
    reader->AddVariable(
    DataVector[ i ].VariableInformation,DataVector[ i ].Data);
    } // END if onlyt metadata
    } // END for all variables in the file

  if( !OnlyMetaData )
  {
  cosmotk::MPIUtilities::Printf(comm,"# Reading Data\n");
  reader->ReadData();
  cosmotk::MPIUtilities::Printf(comm,"# Reading Data [DONE]\n");
  }

  reader->Close();
  delete reader;
  cosmotk::MPIUtilities::Printf(comm,"# Close File [DONE]\n");

  // STEP 3: Loop and print out data according to type
  oss.clear();
  oss.str("");
  oss << std::endl;
  for(int j=0; j < N; ++j )
    {
  if( OnlyMetaData )
    {
    break;
    }

    for(unsigned int i=0; i < DataVector.size(); ++i)
      {
      int type =
        cosmotk::GenericIOUtilities::DetectVariablePrimitiveType(
            DataVector[ i ].VariableInformation );


//      if(type == cosmotk::GENERIC_IO_SHORT_TYPE)
//        {
//        PrintData<short>(DataVector[ i ].Data,j);
//        } // END if short
//      else if( type == cosmotk::GENERIC_IO_LONG_TYPE)
//        {
//        PrintData<long>(DataVector[ i ].Data,j);
//        } // END if long
//      else if( type == cosmotk::GENERIC_IO_LONG_LONG_TYPE)
//        {
//        PrintData<long long>(DataVector[ i ].Data,j);
//        } // END if long long
//      else if( type == cosmotk::GENERIC_IO_INT32_TYPE )
      if( type == cosmotk::GENERIC_IO_INT32_TYPE )
        {
        oss << PrintData<int32_t>(DataVector[ i ].Data,j);
        } // END if int32_t
      else if( type == cosmotk::GENERIC_IO_INT64_TYPE )
        {
        oss << PrintData<int64_t>(DataVector[ i ].Data,j);
        } // END if int64_t
      else if( type == cosmotk::GENERIC_IO_UINT32_TYPE )
        {
        oss << PrintData<uint32_t>(DataVector[ i ].Data,j);
        } // END if uint32_t
      else if( type == cosmotk::GENERIC_IO_UINT64_TYPE )
        {
        oss << PrintData<uint64_t>(DataVector[ i ].Data,j);
        } // END if uint64_t
      else if( type == cosmotk::GENERIC_IO_DOUBLE_TYPE )
        {
        oss << PrintData<double>(DataVector[ i ].Data,j);
        } // END if double
      else if( type == cosmotk::GENERIC_IO_FLOAT_TYPE )
        {
        oss << PrintData<float>(DataVector[ i ].Data,j);
        } // END if float
      else if( type == cosmotk::GENERIC_IO_UCHAR_TYPE )
        {
        oss << PrintData<unsigned char>(DataVector[ i ].Data, j);
        }
      else
        {
        std::cerr << "ERROR: Undefined datatype!\n";
        } // END else
      } // END for all variables read
    oss << std::endl;
    } // END for all elements

  cosmotk::MPIUtilities::SynchronizedPrintf(comm,"%s",oss.str().c_str());
  MPI_Barrier(comm);

  // STEP 4: Delete all read data
  for(unsigned int i=0; i < DataVector.size(); ++i )
    {
    delete [] static_cast<char*>(DataVector[i].Data);
    }

  // STEP 5: Finalize
  MPI_Finalize();
  return 0;
}


