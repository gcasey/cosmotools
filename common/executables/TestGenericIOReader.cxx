/**
 * @brief Simple code to test the GenericIOMPIReader
 */

#include "CosmologyToolsMacros.h"
#include "GenericIOMPIReader.h"

// C/C++ includes
#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <iomanip>

int main(int argc, char **argv)
{

  MPI_Init(&argc,&argv);

  cosmotk::GenericIOMPIReader *myReader =
    new cosmotk::GenericIOMPIReader();
  myReader->SetCommunicator(MPI_COMM_WORLD);

  std::string fileName = std::string(argv[1]);
  myReader->SetFileName(fileName);
  myReader->OpenAndReadHeader();

  std::cout << "Total Number of Blocks: "
       << myReader->GetTotalNumberOfBlocks() << std::endl;
  std::cout << "Assigned blocks: "
       << myReader->GetNumberOfAssignedBlocks() << std::endl;
  std::cout << "Total Number of Elements: "
       << myReader->GetTotalNumberOfElements() << std::endl;
  std::cout.flush();

  std::cout << "Calling GetNumber of Elements....\n";
  std::cout.flush();
  int NumElements = myReader->GetNumberOfElements();
  std::cout << "NumElements: " << NumElements << std::endl;
  std::cout.flush();
  MPI_Barrier(MPI_COMM_WORLD);

  ID_T *haloID = new ID_T[NumElements];
  POSVEL_T *center_x = new POSVEL_T[NumElements];
  POSVEL_T *center_y = new POSVEL_T[NumElements];
  POSVEL_T *center_z = new POSVEL_T[NumElements];
  POSVEL_T *halo_vx  = new POSVEL_T[NumElements];
  POSVEL_T *halo_vy  = new POSVEL_T[NumElements];
  POSVEL_T *halo_vz  = new POSVEL_T[NumElements];
  POSVEL_T *halomass = new POSVEL_T[NumElements];

  myReader->AddVariable("fof_halo_tag",haloID);
  myReader->AddVariable("fof_halo_center_x",center_x);
  myReader->AddVariable("fof_halo_center_y",center_y);
  myReader->AddVariable("fof_halo_center_z",center_z);
  myReader->AddVariable("fof_halo_mean_vx",halo_vx);
  myReader->AddVariable("fof_halo_mean_vy",halo_vy);
  myReader->AddVariable("fof_halo_mean_vz",halo_vz);
  myReader->AddVariable("fof_halo_mass",halomass);

  myReader->ReadData();
  MPI_Barrier(MPI_COMM_WORLD);

  for( int i=0; i < NumElements; ++i )
    {
    std::cout << haloID[ i ] << "\t";
    std::cout << std::scientific
              << std::setprecision(std::numeric_limits<POSVEL_T>::digits10);
    std::cout << center_x[i] << ", " << center_y[i] << ", " << center_z[i];
    std::cout << "\t";
    std::cout << halo_vx[i] << ", " << halo_vy[i] << ", " << halo_vz[i];
    std::cout << "\t";
    std::cout << halomass[i] << std::endl;
    std::cout.flush();
    } // END for all elements

  delete [] haloID;
  delete [] center_x;
  delete [] center_y;
  delete [] center_z;
  delete [] halo_vx;
  delete [] halo_vy;
  delete [] halo_vz;


  myReader->Close();
  MPI_Finalize();
  return 0;
}



