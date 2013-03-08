// C++ includes
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>

#include "CosmologyToolsMacros.h"

// MPI includes
#include <mpi.h>

#ifdef USEDAX
#define DAX_DEVICE_ADAPTER DAX_DEVICE_ADAPTER_TBB
#include "dax/cont/DeviceAdapter.h"
#include "dax/Types.h"
#include "dax/cont/Scheduler.h"
#include "dax/cont/ArrayHandle.h"
#include "dax/mapPointsToGrid.worklet"
#endif

int main( int argc, char **argv )
{
#ifdef USEDAX
  int N = 10;
  REAL *pnt = new REAL[ 3*N ];

  dax::Vector3 origin(0.0);
  dax::Vector3 spacing(0.5);
  dax::Extent3 extent;
  extent.Min = dax::Id3(0);
  extent.Max = dax::Id3(10);


  dax::cont::ArrayHandle<dax::Vector3> pntHandle =
      dax::cont::make_ArrayHandle(reinterpret_cast<dax::Vector3*>(pnt),N);

  dax::cont::ArrayHandle<dax::Id > output;
  dax::cont::Scheduler<> scheduler;
  scheduler.Invoke(
      dax::worklet::MapPointToGrid( ), pntHandle, origin, spacing, extent, output );
#endif
}

// CosmologyTools includes
//#include "CosmologyToolsMacros.h"
//#include "ExtentPartitioner.h"
//#include "ParallelStructureFormationProbe.h"
//#include "StructureFormationProbe.h"
//
//#define IMIN(ext) ext[0]
//#define IMAX(ext) ext[1]
//#define JMIN(ext) ext[2]
//#define JMAX(ext) ext[3]
//#define KMIN(ext) ext[4]
//#define KMAX(ext) ext[5]
//
////-----------------------------------------------------------------------------
//void GetParticlesFromFile(
//    std::string file, REAL *&particles, int &N )
//{
//  std::ifstream ifs;
//  ifs.open( file.c_str() );
//  assert("pre: cannot open file" && ifs.is_open() );
//
//  int NumParticles = 0;
//  ifs >> NumParticles;
//  std::cout << "NumParticles=" << NumParticles << std:: endl;
//  std::cout.flush();
//
//  if( particles != NULL )
//    {
//    delete [] particles;
//    }
//  particles = new REAL[ 3*NumParticles ];
//  assert("pre: particles is NULL" && (particles != NULL) );
//  for( int particleIdx=0; particleIdx < NumParticles; ++particleIdx )
//    {
//    ifs >> particles[ particleIdx*3   ];
//    ifs >> particles[ particleIdx*3+1 ];
//    ifs >> particles[ particleIdx*3+2 ];
//    } // END for all particles
//  ifs.close();
//}
//
////-----------------------------------------------------------------------------
//void GetParticles(
//    REAL *&particles, int &N, int ext[6], REAL origin[3], REAL h[3])
//{
//  N = 0;
//  int idim = IMAX(ext)-IMIN(ext)+1;
//  int jdim = JMAX(ext)-JMIN(ext)+1;
//  int kdim = KMAX(ext)-KMIN(ext)+1;
//  N        = idim * jdim * kdim;
//  int N2   = jdim;
//  int N1   = idim;
//
//  particles = new REAL[ N*3 ];
//  assert("post: particles should not be NULL!");
//
//  int ijk[3];
//  for( int i=IMIN(ext); i <= IMAX(ext); ++i )
//    {
//    for( int j=JMIN(ext); j <= JMAX(ext); ++j )
//      {
//      for( int k=KMIN(ext); k <= KMAX(ext); ++k )
//        {
//        // Compute local IJK
//        ijk[0] = i-IMIN(ext);
//        ijk[1] = j-JMIN(ext);
//        ijk[2] = k-KMIN(ext);
//
//        // Compute linear index into particles array
//        int idx = (ijk[2]*N2+ijk[1])*N1+ijk[0];
//
//        for(int dim=0; dim < 3; ++dim )
//          {
//          particles[idx*3+dim]=origin[dim]+ijk[dim]*h[dim];
//          }
//
//        } // END for all k
//      } // END for all j
//    } // END for all i
//
//}

//------------------------------------------------------------------------------
//int main(int argc, char **argv)
//{
//  int rank, size;
//  MPI_Init(&argc, &argv);
//  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
//  MPI_Comm_size( MPI_COMM_WORLD, &size );
//
//  int ext[6];
//  IMIN(ext) = atoi(argv[1]);
//  IMAX(ext) = atoi(argv[2]);
//  JMIN(ext) = atoi(argv[3]);
//  JMAX(ext) = atoi(argv[4]);
//  KMIN(ext) = atoi(argv[5]);
//  KMAX(ext) = atoi(argv[6]);
//
//  int myext[6];
//  cosmologytools::ExtentPartitioner partitioner;
//  partitioner.SetGlobalExtent( ext );
//  partitioner.SetNumberOfPartitions( size );
//  partitioner.Partition();
//  partitioner.GetExtent( rank, myext );
//
//  int N = 0;
//  REAL h[3];
//  h[0] = h[1] = h[2] = 1.0;
//  REAL origin[3];
//  for( int i=0; i < 3; ++i )
//    {
//    origin[ i ] = 0.0 + myext[ i*2 ]*h[i];
//    }
//  REAL *particles    = NULL;
//  INTEGER *globalIds = NULL;
//  GetParticles(particles, N, myext,origin,h);
//  MPI_Barrier( MPI_COMM_WORLD );
//
//  cosmologytools::ParallelStructureFormationProbe *mySFProbe =
//      new cosmologytools::ParallelStructureFormationProbe();
//  mySFProbe->SetCommunicator( MPI_COMM_WORLD );
//  mySFProbe->SetWholeExtent( ext );
//  mySFProbe->SetParticles( particles, globalIds, N );
//  mySFProbe->Tesselate();
//  MPI_Barrier( MPI_COMM_WORLD );
//
//  delete [] particles;
//  delete mySFProbe;
//  MPI_Finalize();
//  return 0;
//}




