#include <iostream>
#include <mpi.h>
#include <cassert>
#include <cstdlib>

#include "HACCToolsMacros.h"
#include "ExtentPartitioner.h"
#include "ParallelStructureFormationProbe.h"
#include "StructureFormationProbe.h"

#define IMIN(ext) ext[0]
#define IMAX(ext) ext[1]
#define JMIN(ext) ext[2]
#define JMAX(ext) ext[3]
#define KMIN(ext) ext[4]
#define KMAX(ext) ext[5]

//-----------------------------------------------------------------------------
void GetParticles(
    REAL *&particles, int &N, int ext[6], REAL origin[3], REAL h[3])
{
  N = 0;
  int idim = IMAX(ext)-IMIN(ext)+1;
  int jdim = JMAX(ext)-JMIN(ext)+1;
  int kdim = KMAX(ext)-KMIN(ext)+1;
  N        = idim * jdim * kdim;
  int N2   = jdim;
  int N1   = idim;

  particles = new REAL[ N*3 ];
  assert("post: particles should not be NULL!");

  int ijk[3];
  for( int i=IMIN(ext); i <= IMAX(ext); ++i )
    {
    for( int j=JMIN(ext); j <= JMAX(ext); ++j )
      {
      for( int k=KMIN(ext); k <= KMAX(ext); ++k )
        {
        // Compute local IJK
        ijk[0] = i-IMIN(ext);
        ijk[1] = j-JMIN(ext);
        ijk[2] = k-KMIN(ext);

        // Compute linear index into particles array
        int idx = (ijk[2]*N2+ijk[1])*N1+ijk[0];

        for(int dim=0; dim < 3; ++dim )
          {
          particles[idx*3+dim]=origin[dim]+ijk[dim]*h[dim];
          }

        } // END for all k
      } // END for all j
    } // END for all i

}

//------------------------------------------------------------------------------
int main(int argc, char **argv)
{
  int rank, size;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );

  int ext[6];
  IMIN(ext) = atoi(argv[1]);
  IMAX(ext) = atoi(argv[2]);
  JMIN(ext) = atoi(argv[3]);
  JMAX(ext) = atoi(argv[4]);
  KMIN(ext) = atoi(argv[5]);
  KMAX(ext) = atoi(argv[6]);

  int myext[6];
  hacctools::ExtentPartitioner partitioner;
  partitioner.SetGlobalExtent( ext );
  partitioner.SetNumberOfPartitions( size );
  partitioner.Partition();
  partitioner.GetExtent( rank, myext );

  int N = 0;
  REAL h[3];
  h[0] = h[1] = h[2] = 1.0;
  REAL origin[3];
  for( int i=0; i < 3; ++i )
    {
    origin[ i ] = 0.0 + myext[ i*2 ]*h[i];
    }
  REAL *particles = NULL;
  GetParticles(particles, N, myext,origin,h);
  MPI_Barrier( MPI_COMM_WORLD );

  hacctools::ParallelStructureFormationProbe *mySFProbe =
      new hacctools::ParallelStructureFormationProbe();
  mySFProbe->SetCommunicator( MPI_COMM_WORLD );
  mySFProbe->SetWholeExtent( ext );
  mySFProbe->SetParticles( particles );
  mySFProbe->SetNumParticles( N );
  mySFProbe->SetStride( 3 );
  mySFProbe->Tesselate();
  MPI_Barrier( MPI_COMM_WORLD );

  delete [] particles;
  delete mySFProbe;
  MPI_Finalize();
  return 0;
}




