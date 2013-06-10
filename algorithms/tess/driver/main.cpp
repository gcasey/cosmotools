#include "mpi.h"
#include <assert.h>
#include "tess.h"

void GetArgs(int argc, char **argv, int &tb, int *dsize, float *jitter,
	     float *cell_size, float *ghost_factor, float *minvol,
	     float *maxvol);

int main(int argc, char *argv[]) {

  int tb; // total number of blocks in the domain
  int dsize[3]; // domain grid size
  float jitter; // max amount to randomly displace particles
  float cell_size; // max expected cell diameter
  float ghost_factor; // ghost size multiplier
  float minvol, maxvol; // volume range, -1.0 = unused
  double times[MAX_TIMES]; // timing

  MPI_Init(&argc, &argv);

  GetArgs(argc, argv, tb, dsize, &jitter, &cell_size, &ghost_factor,
	  &minvol, &maxvol);
  tess_test(tb, dsize, jitter, cell_size, ghost_factor, minvol, maxvol, times);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Finalize();

  return 0;

}
//----------------------------------------------------------------------------
//
// gets command line args
//
void GetArgs(int argc, char **argv, int &tb, int *dsize, float *jitter,
	     float *cell_size, float *ghost_factor, float *minvol,
	     float *maxvol) {

  assert(argc >= 10);

  tb = atoi(argv[1]);
  dsize[0] = atoi(argv[2]);
  dsize[1] = atoi(argv[3]);
  dsize[2] = atoi(argv[4]);
  *jitter = atof(argv[5]);
  *cell_size = atof(argv[6]);
  *ghost_factor = atof(argv[7]);
  *minvol = atof(argv[8]);
  *maxvol = atof(argv[9]);

}
//----------------------------------------------------------------------------
