/**
 * @brief A simple code to illustrate how to use the GenericIO reader to
 * read in HACC files.
 */

// C/C++ includes
#include <cassert>
#include <iostream>
#include <vector>

// Include MPI
#include <mpi.h>

// CosmoTools includes
#include "CosmologyToolsDefinitions.h"
#include "GenericIODefinitions.hpp"
#include "GenericIOMPIReader.h"
#include "GenericIOReader.h"
#include "GenericIOUtilities.h"

// Global variables
MPI_Comm comm = MPI_COMM_WORLD;
int Rank;
int NumRanks;
int NumElements; /* Number of elements to read in this process */

struct GenericIOData
{
  int GenericIOType;
  cosmotk::VariableInfo VariableInformation;
  void *Data;
};

//==============================================================================
// Global Methods
//==============================================================================
template<class T>
void ProcessData( void *dataPtr, const int NumElements )
{
  assert("pre: input data is NULL!" && (dataPtr != NULL) );
  T* data = static_cast< T* >( dataPtr );

  // TODO: Add your processing here
  for(int i=0; i < NumElements; ++i )
    {
    std::cout << data[ i ] << std::endl;
    } // END loop through all elements
}

/**
 * @brief Program main
 * @param argc argument counter
 * @param argv argument vector
 * @return rc return code
 */
int main(int argc, char** argv)
{
  // STEP 0: Initialize MPI
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(comm,&Rank);
  MPI_Comm_size(comm,&NumRanks);

  // STEP 1: Get file name to read from cmd
  std::string fileName = std::string(argv[1]);

  // STEP 2: Initialize the reader
  cosmotk::GenericIOMPIReader *particleReader =
      new cosmotk::GenericIOMPIReader();
  particleReader->SetCommunicator(comm);
  particleReader->SetFileName(fileName);

  // STEP 3: Open and read header
  particleReader->OpenAndReadHeader();

  // STEP 4: Get the number of elements to read in this process
  NumElements = particleReader->GetNumberOfElements();

  // STEP 5: Register to read all variables from the file
  std::vector< GenericIOData > dataVector;
  dataVector.resize( particleReader->GetNumberOfVariablesInFile() );
  for(int var=0; var < particleReader->GetNumberOfVariablesInFile(); ++var)
    {
    // Get meta-data of the variable
    dataVector[ var ].VariableInformation =
        particleReader->GetFileVariableInfo( var );

    // Get the corresponding primitive type
    dataVector[ var ].GenericIOType =
        cosmotk::GenericIOUtilities::DetectVariablePrimitiveType(
            dataVector[ var ].VariableInformation );
    assert("pre: undefined data type!" &&
           (dataVector[ var ].GenericIOType >= 0) &&
           (dataVector[ var ].GenericIOType < cosmotk::NUM_PRIMITIVE_TYPES) );

    // Allocate buffer to read in the data
    dataVector[ var ].Data =
        cosmotk::GenericIOUtilities::AllocateVariableArray(
            dataVector[ var ].VariableInformation,NumElements);

    // Register variable to be read in the given buffer
    particleReader->AddVariable(
        dataVector[ var ].VariableInformation, dataVector[ var ].Data );
    } // END for all variables


  // STEP 6: Read all registered variables
  particleReader->ReadData();

  // STEP 7: Close the reader
  particleReader->Close();
  delete particleReader;

  // STEP 8: Process the data
  for(unsigned int i=0; i < dataVector.size(); ++i)
    {
    if( dataVector[i].GenericIOType == cosmotk::GENERIC_IO_INT32_TYPE)
      {
      ProcessData<int32_t>(dataVector[i].Data,NumElements);
      } // END if
    else if( dataVector[i].GenericIOType == cosmotk::GENERIC_IO_INT64_TYPE)
      {
      ProcessData<int64_t>(dataVector[i].Data,NumElements);
      } // END else if
    else if( dataVector[i].GenericIOType == cosmotk::GENERIC_IO_UINT32_TYPE )
      {
      ProcessData<uint32_t>(dataVector[i].Data,NumElements);
      } // END else if
    else if( dataVector[i].GenericIOType == cosmotk::GENERIC_IO_UINT64_TYPE )
      {
      ProcessData<uint64_t>(dataVector[i].Data,NumElements);
      } // END else if
    else if( dataVector[i].GenericIOType == cosmotk::GENERIC_IO_DOUBLE_TYPE )
      {
      ProcessData<double>(dataVector[i].Data,NumElements);
      } // END else if
    else if(dataVector[i].GenericIOType == cosmotk::GENERIC_IO_FLOAT_TYPE )
      {
      ProcessData<float>(dataVector[i].Data,NumElements);
      } // END else if
    else
      {
      std::cerr << "ERROR: Undefined datatype!\n";
      MPI_Abort(comm,-1);
      } // END else
    } // END for all variables

  // STEP 9: Clean up all dynamically allocated memory
  for( unsigned int i=0; i < dataVector.size(); ++i)
    {
    delete [] static_cast<char*>(dataVector[i].Data);
    } // END for all variables

  // STEP 10: Finalize MPI
  MPI_Finalize();
  return 0;
}
