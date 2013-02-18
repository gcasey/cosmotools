/**
 * @brief Simple code to test the GenericIOMPIReader
 */

#include <iostream>

#include "GenericIO.h"

int Rank;
int NElem;

int main(int argc, char **argv)
{
  MPI_Init(&argc,&argv);
  MPI_Comm_rank(MPI_COMM_WORLD,&Rank);

  NElem = 10 + Rank*10;
  double *data = new double[NElem];
  int *idx = new int[NElem];
  for( int i=0; i < NElem; ++i)
    {
    data[ i ] = Rank*i + 0.5;
    idx[ i ] = Rank;
    } // END for all elements

  cosmotk::GenericIO writer(MPI_COMM_WORLD,"fake.dat");
  writer.setNumElems(NElem);
  writer.addVariable("DATA",data);
  writer.addVariable("ID",idx);
  writer.write();

  delete [] data;
  delete [] idx;

  MPI_Finalize();
  return 0;
}






