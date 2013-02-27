/**
 * @brief Simple code to test the GenericIOMPIReader
 */

#include <iostream>

#include "GenericIOMPIReader.h"

int main(int argc, char **argv)
{

  MPI_Init(&argc,&argv);

  cosmotk::GenericIOMPIReader *myReader =
    new cosmotk::GenericIOMPIReader();
  myReader->SetCommunicator(MPI_COMM_WORLD);
  myReader->SetFileName( std::string(argv[1]) );

  myReader->OpenAndReadHeader();
  std::cout << "Total Number of Blocks: "
            << myReader->GetTotalNumberOfBlocks() << std::endl;
  std::cout << "Assigned blocks: "
            << myReader->GetNumberOfAssignedBlocks() << std::endl;
  std::cout << "Total Number of Elements: "
            << myReader->GetTotalNumberOfElements() << std::endl;
  std::cout.flush();

  int NumElements = myReader->GetNumberOfElements();
  std::cout << "NumElements: " << NumElements << std::endl;
  std::cout.flush();

  for(int i=0; i < myReader->GetNumberOfVariables(); ++i )
    {
    std::cout << myReader->GetVariableName( i ) << std::endl;
    std::cout.flush();
    }
//  myReader->ReadData();

//  std::cout << "Read: " << NumElements << std::endl;
//  for( int i=0; i < NumElements; ++i)
//    {
//    std::cout  << idx[i] << "\t" << data[i] << std::endl;
//    std::cout.flush();
//    } // END for all elements

  myReader->Close();
  MPI_Finalize();
  return 0;
}



