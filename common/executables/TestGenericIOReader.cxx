/**
 * @brief Simple code to test the GenericIOMPIReader
 */

// C/C++ includes
#include <iostream>
#include <map>
#include <set>
#include <sstream>


#include "GenericIOMPIReader.h"

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

  int NumElements = myReader->GetNumberOfElements();
  std::cout << "NumElements: " << NumElements << std::endl;
  std::cout.flush();

  for(int i=0; i < myReader->GetNumberOfVariables(); ++i )
    {
    std::cout << myReader->GetVariableName( i ) << std::endl;
    std::cout.flush();
    }

  int *rank = new int[ NumElements ];
  int *part = new int[ NumElements ];
  myReader->AddVariable("$rank",rank);
  myReader->AddVariable("$partition",part);
  myReader->ReadData();

  std::map< int,int > fileMap;
  std::set< int > fileSuffix;
  std::cout << "Read: " << NumElements << std::endl;
  std::cout.flush();
  for( int i=0; i < NumElements; ++i)
    {
    fileMap[ rank[i] ] = part[i];
    fileSuffix.insert( part[i] );
    } // END for all elements

  std::cout << "Num files: " << fileSuffix.size() << std::endl;
  std::cout.flush();

  std::ostringstream oss;
  std::vector< cosmotk::GenericIOMPIReader*  > Readers;
  Readers.resize( fileSuffix.size() );
  for(int i=0; i < fileSuffix.size(); ++i)
    {
    oss.str(""); oss.clear();
    oss << fileName << "#" << i;
    std::cout << "========================\n";
    std::cout << "Reading file: " << oss.str() << std::endl;
    std::cout.flush();

    Readers[ i ] = new cosmotk::GenericIOMPIReader();
    Readers[ i ]->SetCommunicator(MPI_COMM_WORLD);
    Readers[ i ]->SetFileName( oss.str() );
    Readers[ i ]->OpenAndReadHeader();

    std::cout << "Total Number of Blocks: "
         << Readers[ i ]->GetTotalNumberOfBlocks() << std::endl;
    std::cout << "Assigned blocks: "
         << Readers[ i ]->GetNumberOfAssignedBlocks() << std::endl;
    std::cout << "Total Number of Elements: "
         << Readers[ i ]->GetTotalNumberOfElements() << std::endl;
    std::cout.flush();

    int NumElements = Readers[ i ]->GetNumberOfElements();
    std::cout << "NumElements: " << NumElements << std::endl;
    std::cout.flush();

    for(int j=0; j < Readers[ i ]->GetNumberOfVariables(); ++j )
      {
      std::cout << Readers[ i ]->GetVariableName( j ) << std::endl;
      std::cout.flush();
      }

    Readers[ i ]->Close();
    delete Readers[ i ];
    } // END for all files

  delete [] rank;
  delete [] part;


  myReader->Close();
  MPI_Finalize();
  return 0;
}



