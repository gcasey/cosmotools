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

  myReader->Close();

  MPI_Finalize();
  return 0;
}



