/**
 * @brief A simple program that opens a file written in the GenericIO format
 * and prints out the number of variables, their corresponding sizes and
 * names.
 */

// C++ includes
#include <iostream>
#include <vector>
#include <iomanip>
#include <limits>

// MPI include
#include <mpi.h>

// CosmologyTools includes
#include "GenericIODefinitions.hpp"
#include "GenericIOMPIReader.h"
#include "GenericIOUtilities.h"
#include "MPIUtilities.h"

//==============================================================================
// Global variables
//==============================================================================
int rank;
int size;
MPI_Comm comm = MPI_COMM_WORLD;

struct DataInformation
{
  cosmotk::VariableInfo VariableInformation;
  void *Data;
};

//==============================================================================
// Global Methods
//==============================================================================
template<class T>
void PrintData(void *dataPtr, const int idx)
{
  T* castedPtr = static_cast< T* >( dataPtr );
  std::cout << std::scientific
            << std::setprecision( std::numeric_limits<T>::digits10 )
            << castedPtr[ idx ] << ";";
  std::cout.flush();
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
//  cosmotk::MPIUtilities::Printf(comm,"- Initialize MPI...[DONE]\n");

  // STEP 1: Get file name to read
  std::string file=std::string(argv[1]);

  // STEP 2: Read and print vars
  cosmotk::GenericIOMPIReader reader;
  reader.SetCommunicator(comm);
  reader.SetFileName(file);
  reader.OpenAndReadHeader();

//  cosmotk::MPIUtilities::SynchronizedPrintf(
//      comm,"NumElements: %d\n", reader.GetNumberOfElements() );
//  cosmotk::MPIUtilities::SynchronizedPrintf(
//      comm,"NumVariables: %d\n", reader.GetNumberOfVariablesInFile() );


  std::vector< DataInformation > DataVector;
  DataVector.resize( reader.GetNumberOfVariablesInFile() );
  int N = reader.GetNumberOfElements();
  for(int i=0; i < reader.GetNumberOfVariablesInFile(); ++i)
    {
//    cosmotk::MPIUtilities::Printf(
//        comm,"var[%d]=%s\n", i,
//        const_cast<char*>(reader.GetVariableName(i).c_str()) );
    std::cout << reader.GetVariableName(i) << ";";
    DataVector[ i ].VariableInformation = reader.GetFileVariableInfo( i );
    DataVector[ i ].Data =
        cosmotk::GenericIOUtilities::AllocateVariableArray(
            DataVector[ i ].VariableInformation,N);
    reader.AddVariable(
        DataVector[ i ].VariableInformation,DataVector[ i ].Data);
    }
  std::cout << std::endl;
  reader.ReadData();
  reader.Close();

  // STEP 3: Loop and print out data according to type
  for(int j=0; j < N; ++j )
    {
    for(unsigned int i=0; i < DataVector.size(); ++i)
      {
      int type =
        cosmotk::GenericIOUtilities::DetectVariablePrimitiveType(
            DataVector[ i ].VariableInformation );


      if(type == cosmotk::GENERIC_IO_SHORT_TYPE)
        {
        PrintData<short>(DataVector[ i ].Data,j);
        } // END if short
      else if( type == cosmotk::GENERIC_IO_LONG_TYPE)
        {
        PrintData<long>(DataVector[ i ].Data,j);
        } // END if long
      else if( type == cosmotk::GENERIC_IO_LONG_LONG_TYPE)
        {
        PrintData<long long>(DataVector[ i ].Data,j);
        } // END if long long
      else if( type == cosmotk::GENERIC_IO_INT32_TYPE )
        {
        PrintData<int32_t>(DataVector[ i ].Data,j);
        } // END if int32_t
      else if( type == cosmotk::GENERIC_IO_INT64_TYPE )
        {
        PrintData<int64_t>(DataVector[ i ].Data,j);
        } // END if int64_t
      else if( type == cosmotk::GENERIC_IO_UINT32_TYPE )
        {
        PrintData<uint32_t>(DataVector[ i ].Data,j);
        } // END if uint32_t
      else if( type == cosmotk::GENERIC_IO_UINT64_TYPE )
        {
        PrintData<uint64_t>(DataVector[ i ].Data,j);
        } // END if uint64_t
      else if( type == cosmotk::GENERIC_IO_DOUBLE_TYPE )
        {
        PrintData<double>(DataVector[ i ].Data,j);
        } // END if double
      else if( type == cosmotk::GENERIC_IO_FLOAT_TYPE )
        {
        PrintData<float>(DataVector[ i ].Data,j);
        } // END if float
      else
        {
        std::cerr << "ERROR: Undefined datatype!\n";
        } // END else
      } // END for all variables read
    std::cout << std::endl;
    } // END for all elements

  // STEP 4: Delete all read data
  for(unsigned int i=0; i < DataVector.size(); ++i )
    {
    delete [] static_cast<char*>(DataVector[i].Data);
    }

  // STEP 5: Finalize
  MPI_Finalize();
  return 0;
}


