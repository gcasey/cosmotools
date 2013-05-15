/*---------------------------------------------------------------------------
 *
 * parallel voronoi tesselation
 *
 * Tom Peterka
 * Argonne National Laboratory
 * 9700 S. Cass Ave.
 * Argonne, IL 60439
 * tpeterka@mcs.anl.gov
 *
 * (C) 2013 by Argonne National Laboratory.
 * See COPYRIGHT in top-level directory.
 *
--------------------------------------------------------------------------*/
#include "mpi.h"
#include "diy.h"
#include "tess.h"
#include "io.h"

#include <stddef.h>
#include <stdio.h>

static int dim = 3; /* everything 3D */
static float data_mins[3], data_maxs[3]; /* extents of overall domain */
static MPI_Comm comm; /* MPI communicator */
static float min_vol, max_vol; /* cell volume range */

static int nblocks; /* number of blocks per process */
static float cell_size; /* maximum expected cell diameter */
static float ghost_factor; /* ghost size multiplier */
static double *times; /* timing info */

/* datatype profiling for John Jenkins' research */
/* #define DATATYPE_PROFILE */

#ifdef DATATYPE_PROFILE
static int **hdrs; /* headers */
#endif

/* #define PNETCDF_IO */

/*------------------------------------------------------------------------*/
/*
  initialize parallel voronoi (and later possibly delaunay) tesselation

  num_blocks: local number of blocks in my process
  gids: global ids of my local blocks
  bounds: block bounds (extents) of my local blocks
  neighbors: neighbor lists for each of my local blocks, in lid order
  neighbor bounds need not be known, will be discovered automatically
  num_neighbors: number of neighbors for each of my local blocks, in lid order
  cell_dia: maximum expeted cell diameter
  ghost_mult: ghost size multiplier
  global_mins, global_maxs: overall data extents
  wrap: whether wraparound neighbors are used
  minvol, maxvol: filter range for which cells to keep
  pass -1.0 to skip either or both bounds
  mpi_comm: MPI communicator
  all_times: times for particle exchange, voronoi cells, convex hulls, and output
*/
void tess_init(int num_blocks, int *gids, 
	       struct bb_t *bounds, struct gb_t **neighbors, 
	       int *num_neighbors, float cell_dia, float ghost_mult, 
	       float *global_mins, float *global_maxs, 
	       int wrap, float minvol, float maxvol, MPI_Comm mpi_comm,
	       double *all_times) {

  int i;

  /* save globals */
  comm = mpi_comm;
  nblocks = num_blocks;
  cell_size = cell_dia;
  ghost_factor = ghost_mult;
  min_vol = minvol;
  max_vol = maxvol;
  times = all_times;

  /* data extents */
  for(i = 0; i < 3; i++) {
    data_mins[i] = global_mins[i];
    data_maxs[i] = global_maxs[i];
  }

  /* init times */
  for (i = 0; i < MAX_TIMES; i++)
    times[i] = 0.0;

  /* init DIY */
  DIY_Init(dim, NULL, 1, comm);
  DIY_Decomposed(num_blocks, gids, bounds, NULL, NULL, NULL, NULL, neighbors, 
		 num_neighbors, wrap);

}
/*------------------------------------------------------------------------*/
/*
finalize parallel voronoi (and later possibly delaunay) tesselation
*/
void tess_finalize() {

  DIY_Finalize();

}
/*------------------------------------------------------------------------*/
/*
  tessellate parallel voronoi (and later possibly delaunay) tesselation

  particles: particles[block_num][particle] 
  where each particle is 3 values, px, py, pz
  num_particles; number of particles in each block
  out_file: output file name
*/
void tess(float **particles, int *num_particles, char *out_file) {

  voronoi(nblocks, particles, num_particles, cell_size, ghost_factor,
	  times, out_file);

}
/*------------------------------------------------------------------------*/
/*
  test of parallel voronoi (and later possibly delaunay) tesselation

  tot_blocks: total number of blocks in the domain
  data_size: domain grid size (x, y, z)
  jitter: maximum amount to randomly move each particle
  cell_size: maximum expeted cell diameter
  ghost_factor: ghost size multiplier
  minvol, maxvol: filter range for which cells to keep
  pass -1.0 to skip either or both bounds
  times: times for particle exchange, voronoi cells, convex hulls, and output
*/
void tess_test(int tot_blocks, int *data_size, float jitter, float cell_size,
	       float ghost_factor, float minvol, float maxvol, double *times) {

  float **particles; /* particles[block_num][particle] 
			 where each particle is 3 values, px, py, pz */
  int *num_particles; /* number of particles in each block */
  int dim = 3; /* 3D */
  int given[3] = {0, 0, 0}; /* no constraints on decomposition in {x, y, z} */
  int ghost[6] = {0, 0, 0, 0, 0, 0}; /* ghost in {-x, +x, -y, +y, -z, +z} */
  int nblocks; /* my local number of blocks */
  int i;

  comm = MPI_COMM_WORLD;
  min_vol = minvol;
  max_vol = maxvol;

  /* data extents */
  for(i = 0; i < 3; i++) {
    data_mins[i] = 0.0;
    data_maxs[i] = data_size[i] - 1.0;
  }

  /* have DIY do the decomposition */
  DIY_Init(dim, data_size, 1, comm);
  DIY_Decompose(ROUND_ROBIN_ORDER, tot_blocks, &nblocks, 1, ghost, given);

  /* generate test points in each block */
  particles = (float **)malloc(nblocks * sizeof(float *));
  num_particles = (int *)malloc(nblocks * sizeof(int));
  for (i = 0; i < nblocks; i++)
    num_particles[i] = gen_particles(i, &particles[i], jitter);

  voronoi(nblocks, particles, num_particles, cell_size, ghost_factor,
	  times, "vor.out");

  /* cleanup */
  for (i = 0; i < nblocks; i++)
    free(particles[i]);
  free(particles);
  free(num_particles);

  DIY_Finalize();

}
/*--------------------------------------------------------------------------*/
/*
  postprocessing parallel voronoi (and later possibly delaunay) tesselation

  num_blocks: local number of blocks in my process
  gids: global ids of my local blocks
  bounds: block bounds (extents) of my local blocks
  neighbors: neighbor lists for each of my local blocks, in lid order
  neighbor bounds need not be known, will be discovered automatically
  num_neighbors: number of neighbors for each of my local blocks, in lid order
  cell_dia: maximum expeted cell diameter
  ghost_mult: ghost size multiplier
  global_mins, global_maxs: overall data extents
  wrap: whether wraparound neighbors are used
  minvol, maxvol: filter range for which cells to keep
  pass -1.0 to skip either or both bounds
  mpi_comm: MPI communicator
  all_times: times for particle exchange, voronoi cells, convex hulls, 
  and output
  particles: particles[block_num][particle] 
  where each particle is 3 values, px, py, pz
  num_particles; number of particles in each block
  out_file: output file name
*/
void tess_post(int num_blocks, int *gids, 
	       struct bb_t *bounds, struct gb_t **neighbors, 
	       int *num_neighbors, float cell_dia, float ghost_mult,
	       float *global_mins, float *global_maxs, 
	       int wrap, float minvol, float maxvol, MPI_Comm mpi_comm,
	       double *all_times, float **particles, int *num_particles, 
	       char *out_file) {

  int dim = 3; // 3D
  int i;

  /* save globals */
  comm = mpi_comm;
  nblocks = num_blocks;
  cell_size = cell_dia;
  ghost_factor = ghost_mult;
  min_vol = minvol;
  max_vol = maxvol;
  times = all_times;

  /* data extents */
  for(i = 0; i < 3; i++) {
    data_mins[i] = global_mins[i];
    data_maxs[i] = global_maxs[i];
  }

  /* init times */
  for (i = 0; i < MAX_TIMES; i++)
    times[i] = 0.0;

  /* init DIY */
  DIY_Init(dim, NULL, 1, comm);
  DIY_Decomposed(num_blocks, gids, bounds, NULL, NULL, NULL, NULL, neighbors, 
		 num_neighbors, wrap);

  voronoi(nblocks, particles, num_particles, cell_dia, ghost_mult,
	  times, out_file);

  DIY_Finalize();

}
/*--------------------------------------------------------------------------*/
/*
  test of parallel voronoi (and later possibly delaunay) tesselation

  nblocks: local number of blocks
  particles: particles[block_num][particle] 
  where each particle is 3 values, px, py, pz
  num_particles; number of particles in each block
  cell_size: maximum expeted cell diameter
  ghost_factor: ghost size multiplier
  times: times for particle exchange, voronoi cells, convex hulls, and output
  out_file: output file name
*/
void voronoi(int nblocks, float **particles, int *num_particles, 
	     float cell_size, float ghost_factor, double *times, 
	     char *out_file) {

  int *num_orig_particles; /* number of original particles, before any
			      neighbor exchange */

  int dim = 3; /* 3D */

#ifndef DATATYPE_PROFILE
  int **hdrs; /* headers */
#endif

  struct vblock_t *vblocks; /* voronoi blocks */
  int i;

  num_orig_particles = (int *)malloc(nblocks * sizeof(int));
  for (i = 0; i < nblocks; i++)
    num_orig_particles[i] = num_particles[i];

#ifdef TIMING
  MPI_Barrier(comm);
  times[EXCH_TIME] = MPI_Wtime();
#endif

  /* exchange particles with neighbors */
  neighbor_particles(nblocks, particles,
		     num_particles, cell_size * ghost_factor);

#ifdef TIMING
  MPI_Barrier(comm);
  times[EXCH_TIME] = MPI_Wtime() - times[EXCH_TIME];
#endif

  /* allocate and initialize blocks */
  create_blocks(nblocks, &vblocks, &hdrs);

#ifdef TIMING
  MPI_Barrier(comm);
  times[CELL_TIME] = MPI_Wtime();
#endif

  /* create original voronoi cells */
  orig_cells(nblocks, vblocks, dim, num_particles, num_orig_particles,
	     particles, cell_size * ghost_factor);

#ifdef TIMING
  /* no barrier here; want min and max time */
  times[CELL_TIME] = MPI_Wtime() - times[CELL_TIME];
  MPI_Barrier(comm);
  times[HULL_TIME] = MPI_Wtime();
#endif

  /* create convex hulls for completed cells */
  cell_hulls(nblocks, vblocks, dim);

#ifdef TIMING
  /* no barrier here; want min and max time */
  times[HULL_TIME] = MPI_Wtime() - times[HULL_TIME];
#endif

  /* prepare for output */
  prep_out(nblocks, vblocks, hdrs);

  /* debug */
/*   int b; */
/*   for (b = 0; b < nblocks; b++) { */
/*     struct vblock_t *v = &vblocks[b]; */
/*     fprintf(stderr, "mem gid = %d num_verts = %d num_complete_cells = %d " */
/* 	    "tot_num_cell_faces = %d tot_num_face_verts = %d " */
/* 	    "num_orig_particles = %d\n", DIY_Gid(0, b), */
/* 	    v->num_verts, v->num_complete_cells, v->tot_num_cell_faces, */
/* 	    v->tot_num_face_verts, v->num_orig_particles); */
/*   } */

#ifdef TIMING
  MPI_Barrier(comm);
  times[OUT_TIME] = MPI_Wtime();
#endif

  /* write output */
#ifdef PNETCDF_IO
  char out_ncfile[256];
  strncpy(out_ncfile, out_file, sizeof(out_ncfile));
  strncat(out_ncfile, ".nc", sizeof(out_file));
  pnetcdf_write(nblocks, vblocks, out_ncfile, comm);
#else
  diy_write(nblocks, vblocks, hdrs, out_file);
#endif

#ifdef TIMING
  MPI_Barrier(comm);
  times[OUT_TIME] = MPI_Wtime() - times[OUT_TIME];
#endif

  /* collect stats */
  collect_stats(nblocks, vblocks, times);

  /* cleanup */
  destroy_blocks(nblocks, vblocks, hdrs);
  free(num_orig_particles);

}
/*--------------------------------------------------------------------------*/
/*
  creates original voronoi cells

  nblocks: number of blocks
  vblocks: pointer to array of vblocks
  dim: number of dimensions (eg. 3)
  num_particles: number of particles in each block
  num_orig_particles: number of original particles in each block, before any
  neighbor exchange
  particles: particles in each block, particles[block_num][particle]
  where each particle is 3 values, px, py, pz
  ghost: extra size per side beyond block bounds to consider vertices as valid
*/
void orig_cells(int nblocks, struct vblock_t *vblocks, int dim,
		int *num_particles, int *num_orig_particles, 
		float **particles, float ghost) {

  boolT ismalloc = False;    /* True if qhull should free points in
				qh_freeqhull() or reallocation */
  char flags[250];          /* option flags for qhull, see qh-quick.htm */
  int exitcode;             /* 0 if no error from qhull */
  int curlong, totlong;     /* memory remaining after qh_memfreeshort */
  FILE *dev_null; /* file descriptor for writing to /dev/null */
  int i, j;

  dev_null = fopen("/dev/null", "w");
  assert(dev_null != NULL);

  /* for all blocks */
  for (i = 0; i < nblocks; i++) {

    /* deep copy from float to double (qhull API is double) */
    double *pts = (double *)malloc(num_particles[i] * 3 * sizeof(double));
    for (j = 0; j < 3 * num_particles[i]; j++)
      pts[j] = particles[i][j];

    /* compute voronoi */
    sprintf (flags, "qhull v o Fv"); /* voronoi cells and adjacencies */

    /* eat qhull output by sending it to dev/null
       need to see how this behaves on BG/P, will get I/O forwarded but will 
       stop there and not proceed to storage */
    exitcode = qh_new_qhull(dim, num_particles[i], pts, ismalloc,
			    flags, dev_null, stderr);

    free(pts);

    /* process voronoi output */
    if (!exitcode)
      gen_voronoi_output(qh facet_list, &vblocks[i], num_particles[i]);

    /* allocate cell sites for original particles */
    vblocks[i].num_orig_particles = num_orig_particles[i];
    vblocks[i].sites =
      (float *)malloc(3 * sizeof(float) * vblocks[i].num_orig_particles);
    for (j = 0; j < vblocks[i].num_orig_particles; j++) {
      vblocks[i].sites[3 * j] = particles[i][3 * j];
      vblocks[i].sites[3 * j + 1] = particles[i][3 * j + 1];
      vblocks[i].sites[3 * j + 2] = particles[i][3 * j + 2];
    }
    vblocks[i].num_cells = vblocks[i].num_orig_particles;

    /* determine complete cells */
    complete_cells(&vblocks[i], i, ghost);

    /* clean up qhull */
    qh_freeqhull(!qh_ALL);                 /* free long memory */
    qh_memfreeshort(&curlong, &totlong);  /* free short memory */
    if (curlong || totlong)
      fprintf (stderr, "qhull internal warning: did not free %d bytes of "
	       "long memory (%d pieces)\n", totlong, curlong);

  } /* for all blocks */

  fclose(dev_null);

}
/*--------------------------------------------------------------------------*/
/*
  creates convex hulls for completed cells

  nblocks: number of blocks
  vblocks: pointer to array of vblocks
  dim: number of dimensions (eg. 3)
*/
void cell_hulls(int nblocks, struct vblock_t *vblocks, int dim) {

  boolT ismalloc = False;    /* True if qhull should free points in
				qh_freeqhull() or reallocation */
  char flags[250];          /* option flags for qhull, see qh-quick.htm */
  int exitcode;             /* 0 if no error from qhull */
  int curlong, totlong;     /* memory remaining after qh_memfreeshort */
  double *vertices; /* convex hull vertices */
  int *vmap; /* map of convex hull vertex ids back to original voronoi ids */
  struct cblock_t cblock; /* convex hull block */
  FILE *dev_null; /* file descriptor for writing to /dev/null */
  int i, j;

  /* debug */
  double t0, prep_time, proc_time, post_time;
  prep_time = 0.0;
  proc_time = 0.0;
  post_time = 0.0;

  dev_null = fopen("/dev/null", "w");
  assert(dev_null != NULL);

  /* for all blocks */
  for (i = 0; i < nblocks; i++) {

    vblocks[i].num_complete_cells = 0;

    /* for all complete cells in the current block */
    for (j = 0; j < vblocks[i].temp_num_complete_cells; j++) { 

      /* debug */
      t0 = MPI_Wtime();

      /* prepare vertices and compute convex hull */
      int num_verts = 
	vblocks[i].num_cell_verts[vblocks[i].temp_complete_cells[j]];
      vertices = (double *)malloc(num_verts * 3 * sizeof(double));
      vmap = (int *)malloc(num_verts * sizeof(int));
      prep_vertices(&vblocks[i], j, vertices, vmap);

      /* debug */
      prep_time += (MPI_Wtime() - t0);
      t0 = MPI_Wtime();

      sprintf (flags, "qhull o FS"); /* convex hull */
      exitcode = qh_new_qhull(dim, num_verts, vertices, ismalloc,
			      flags, dev_null, stderr);

      /* debug */
      proc_time += (MPI_Wtime() - t0);
      t0 = MPI_Wtime();

      /* compute convex hull, check volume against filter range, and
	 store the complete cell if it passes */
      if (gen_convex_output(qh facet_list, &cblock)) {
	vblocks[i].complete_cells[vblocks[i].num_complete_cells] = 
	  vblocks[i].temp_complete_cells[j];
	vblocks[i].num_complete_cells++;
	convex_to_voronoi(&cblock, &vblocks[i], vmap, 
			  vblocks[i].num_complete_cells - 1);
      }

      /* debug */
      post_time += (MPI_Wtime() - t0);

      /* cleanup */
      free(vertices);
      free(vmap);
      qh_freeqhull(!qh_ALL);                 /* free long memory */
      qh_memfreeshort(&curlong, &totlong);  /* free short memory and */
      if (curlong || totlong)
	fprintf (stderr, "qhull internal warning: did not free %d bytes of "
		 "long memory (%d pieces)\n", totlong, curlong);

    } /* for all complete cells in the current block */

    /* debug: print the cells */
/*     int n = 0; */
/*     int m = 0; */
/*     int k, l; */
/*     fprintf(stderr, "\n%d final cells in local block lid = %d:\n\n", */
/* 	    vblocks[i].num_complete_cells, i); */
/*     for (j = 0; j < vblocks[i].num_complete_cells; j++) { */
/*       for (k = 0; k < vblocks[i].num_cell_faces[j]; k++) { */
/* 	for (l = 0; l < vblocks[i].num_face_verts[n]; l++) { */
/* 	  int v = vblocks[i].face_verts[m]; */
/* 	  fprintf(stderr, "%d [%.1lf %.1lf %.1lf] ", v, */
/* 		  vblocks[i].verts[3 * v], */
/* 		  vblocks[i].verts[3 * v + 1], */
/* 		  vblocks[i].verts[3 * v + 2]); */
/* 	  m++; */
/* 	} */
/* 	n++; */
/* 	fprintf(stderr, "\n"); */
/*       } */
/*       fprintf(stderr, "area = %.1lf volume = %.1lf\n\n", vblocks[i].areas[j], */
/* 	      vblocks[i].vols[j]); */
/*     } */

  } /* for all blocks */

  fclose(dev_null);

  /* debug */
/*   int rank; */
/*   MPI_Comm_rank(MPI_COMM_WORLD, &rank); */
/*   fprintf(stderr, "rank = %d temp_num_complete_cells = %d num_complete_cells = %d prep_time = %.3lf proc_time = %.3lf post_time = %.3lf\n", */
/* 	  rank, vblocks[0].temp_num_complete_cells,  */
/* 	  vblocks[0].num_complete_cells, prep_time, proc_time, post_time); */

}
/*--------------------------------------------------------------------------*/
/*
  exchanges particles with neighbors

  nblocks: number of blocks
  particles: old particles before neighbor exchange
  num_particles: number of new particles in each block (input / output)
  ghost: proximity threshold to neighbors, to determine wheter near enough 
  to send the site to the neighbor
*/
void neighbor_particles(int nblocks, float **particles, int *num_particles,
			float ghost) {

  void ***recv_particles; /* pointers to particles in ecah block 
			     that are received from neighbors */
  int *num_recv_particles; /* number of received particles for each block */
  int i, j;
  struct bb_t bb; // current block bounds

  recv_particles = (void ***)malloc(nblocks * sizeof(void **));
  num_recv_particles = (int *)malloc(nblocks * sizeof(int));

  /* for all blocks */
  for (i = 0; i < nblocks; i++) {

    DIY_Block_bounds(0, i, &bb);

    /* for all particles in the current block */
    for (j = 0; j < num_particles[i]; j++) {

      /* only enqueue particles close enough to block boundary */
      if (particles[i][3 * j]   -   bb.min[0] < ghost ||
	  bb.max[0]   -   particles[i][3 * j] < ghost ||
	  particles[i][3 * j + 1] - bb.min[1] < ghost ||
	  bb.max[1] - particles[i][3 * j + 1] < ghost ||
	  particles[i][3 * j + 2] - bb.min[2] < ghost ||
	  bb.max[2] - particles[i][3 * j + 2] < ghost)

	DIY_Enqueue_item_all_near(0, i, (void *)(&(particles[i][3 * j])),
				  NULL, 3 * sizeof(float),
				  &particles[i][3 * j],
				  ghost, &transform_particle);

    }

  } /* for all blocks */

  /* exchange neighbors */
  DIY_Exchange_neighbors(0, recv_particles, num_recv_particles, 1.0, 
			 &item_type);

  /* copy received particles to particles */
  for (i = 0; i < nblocks; i++) {

    if (num_recv_particles[i]) {

      /* grow space */
      particles[i] = 
	(float *)realloc(particles[i], 
			 (num_particles[i] + num_recv_particles[i]) *
			 3 * sizeof(float));

      /* copy received particles */
      for (j = 0; j < num_recv_particles[i]; j++) { 
	particles[i][3 * num_particles[i]] =
	  ((float*)recv_particles[i][j])[0];
	particles[i][3 * num_particles[i] + 1] =
	  ((float*)recv_particles[i][j])[1];
	particles[i][3 * num_particles[i] + 2] =
	  ((float*)recv_particles[i][j])[2];
	num_particles[i]++;
      }

    }

  }

  /* clean up */
  DIY_Flush_neighbors(0, recv_particles, num_recv_particles, &item_type);
  free(num_recv_particles);
  free(recv_particles);

}
/*--------------------------------------------------------------------------*/
/*
 makes DIY datatype for sending / receiving one item
*/
void item_type(DIY_Datatype *dtype) {

  DIY_Create_vector_datatype(3, 1, DIY_FLOAT, dtype);

}
/*--------------------------------------------------------------------------*/
/*
  collects statistics

  nblocks: number of blocks
  vblocks: pointer to array of vblocks
  times: timing info
*/
void collect_stats(int nblocks, struct vblock_t *vblocks, double *times) {

  int i, j, k, m, n, v, f;
  int tot_num_cell_verts = 0; /* number of cell verts in all local blocks */
  int tot_num_face_verts = 0; /* number of face verts in all local blocks */
  int unique_verts[AVG_CELL_FACES * AVG_FACE_VERTS]; /* unique vertices in 
							one cell */
  int num_unique_verts; /* number of unique vertices in one cell*/
  float vol_bin_width; /* width of a volume histogram bin */
  float dense_bin_width; /* width of a density histogram bin */
  float tot_cell_vol = 0.0; /* sum of cell volumes */
  float tot_cell_dense = 0.0; /* sum of cell densities */
  struct stats_t stats; /* local stats */
  static int first_dense = 1; /* first density value saved */
  int rank;

  MPI_Comm_rank(comm, &rank);

  /* timing range */
  stats.min_cell_time = times[CELL_TIME];
  stats.max_cell_time = times[CELL_TIME];
  stats.min_hull_time = times[HULL_TIME];
  stats.max_hull_time = times[HULL_TIME];

  /* --- first pass: find average number of vertices per cell and 
     volume range --- */

  stats.tot_cells = 0;
  stats.tot_faces = 0;
  stats.tot_verts = 0;
  stats.avg_cell_verts = 0.0;
  stats.avg_cell_faces = 0.0;
  stats.avg_face_verts = 0.0;
  stats.avg_cell_vol   = 0.0;
  stats.avg_cell_dense = 0.0;

  for (i = 0; i < nblocks; i++) { /* for all blocks */

    stats.tot_cells += vblocks[i].num_complete_cells;
    stats.tot_verts += vblocks[i].num_verts;

    /* for all complete cells in the current block */
    f = 0;
    v = 0;
    for (j = 0; j < vblocks[i].num_complete_cells; j++) {

      if (vblocks[i].vols[j] == 0.0)
	fprintf(stderr, "found cell with 0.0 volume--this should not happen\n");

      tot_cell_vol += vblocks[i].vols[j];
      float dense = 0.0;
      if (vblocks[i].vols[j] > 0.0) {
	dense = 1.0 / vblocks[i].vols[j];
	tot_cell_dense += dense;
      }

      stats.tot_faces += vblocks[i].num_cell_faces[j];
      num_unique_verts = 0;

      /* volume range */
      if (i == 0 && j == 0) {
	stats.min_cell_vol = vblocks[i].vols[j];
	stats.max_cell_vol = vblocks[i].vols[j];
      }
      else {
	if (vblocks[i].vols[j] < stats.min_cell_vol)
	  stats.min_cell_vol = vblocks[i].vols[j];
	if (vblocks[i].vols[j] > stats.max_cell_vol)
	  stats.max_cell_vol = vblocks[i].vols[j];
      }

      /* density range */
      if (first_dense && vblocks[i].vols[j] > 0.0) {
	stats.min_cell_dense = dense;
	stats.max_cell_dense = stats.min_cell_dense;
	first_dense = 0;
      }
      else if (vblocks[i].vols[j] > 0.0) {
	if (dense < stats.min_cell_dense)
	  stats.min_cell_dense = dense;
	if (dense > stats.max_cell_dense)
	  stats.max_cell_dense = dense;
      }

      /* for all faces in the current cell */
      for (k = 0; k < vblocks[i].num_cell_faces[j]; k++) {

	tot_num_face_verts += vblocks[i].num_face_verts[f++];

	/* for all verts in the current face */
	for (m = 0; m < vblocks[i].num_face_verts[k]; m++) {

	  /* check if we already counted it */
	  for (n = 0; n < num_unique_verts; n++) {
	    if (vblocks[i].face_verts[v] == unique_verts[n])
	      break;
	  }
	  if (n == num_unique_verts) {
	    unique_verts[n] = vblocks[i].face_verts[v];
	    num_unique_verts++;
	    assert(num_unique_verts < AVG_CELL_FACES * AVG_FACE_VERTS - 1);
	  }

	  v++;

	} /* for all verts */


      } /* for all faces */

      tot_num_cell_verts += num_unique_verts;

    } /* for all complete cells */


  } /* for all blocks */


  /* compute local averages */

  /* debug */
/*   fprintf(stderr, "tot_num_cell_verts = %d tot_num_face_verts = %d tot_cell_vol = %.3lf tot_cell_dense = %.3lf tot_cells = %d tot_faces = %d\n",  */
/* 	  tot_num_cell_verts, tot_num_face_verts, tot_cell_vol, tot_cell_dense, */
/* 	  stats.tot_cells, stats.tot_faces); */

  if (stats.tot_cells) { /* don't divide by 0 */
    stats.avg_cell_verts = tot_num_cell_verts / stats.tot_cells;
    stats.avg_cell_faces = stats.tot_faces / stats.tot_cells;
    stats.avg_face_verts = tot_num_face_verts / stats.tot_faces;
    stats.avg_cell_vol = tot_cell_vol / stats.tot_cells;
    stats.avg_cell_dense = tot_cell_dense / stats.tot_cells;
  }

  /* aggregate totals across all procs and compute average on that */
  aggregate_stats(nblocks, vblocks, &stats);

  /* --- print output --- */

  /* global stats */
  vol_bin_width = (stats.max_cell_vol - stats.min_cell_vol) / 
    stats.num_vol_bins;
  dense_bin_width = (stats.max_cell_dense - stats.min_cell_dense) / 
    stats.num_dense_bins;
  if (rank == 0) {
    fprintf(stderr, "----------------- global stats ------------------\n");
    fprintf(stderr, "particle exchange time = %.3lf s\n", times[EXCH_TIME]);
    fprintf(stderr, "[min, max] voronoi cell time = [%.3lf, %.3lf] s\n",
	    stats.min_cell_time, stats.max_cell_time);
    fprintf(stderr, "[min, max] convex hull time = [%.3lf, %.3lf] s\n",
	    stats.min_hull_time, stats.max_hull_time);
    fprintf(stderr, "output time = %.3lf s\n", times[OUT_TIME]);
    fprintf(stderr, "-----\n");
    fprintf(stderr, "total cells found = %d\n", stats.tot_cells);
    fprintf(stderr, "total cell vertices found = %d\n", stats.tot_verts);
    fprintf(stderr, "average number of vertices per cell = %.0lf\n",
	    stats.avg_cell_verts);
    fprintf(stderr, "average number of faces per cell = %.0lf\n",
	    stats.avg_cell_faces);
    fprintf(stderr, "average number of vertices per face = %.0lf\n",
	    stats.avg_face_verts);
    fprintf(stderr, "-----\n");
    fprintf(stderr, "min cell volume = %.3lf max cell volume = %.3lf "
	    "avg cell volume = %.3lf units^3\n",
	    stats.min_cell_vol, stats.max_cell_vol, stats.avg_cell_vol);
    fprintf(stderr, "number of cell volume histogram bins = %d\n",
	    stats.num_vol_bins);
    fprintf(stderr, "-----\n");
    fprintf(stderr, "cell volume histogram:\n");
    fprintf(stderr, "min value\tcount\t\tmax value\n");
    for (k = 0; k < stats.num_vol_bins; k++)
      fprintf(stderr, "%.3lf\t\t%d\t\t%.3lf\n", 
	      stats.min_cell_vol + k * vol_bin_width, stats.vol_hist[k], 
	      stats.min_cell_vol + (k + 1) * vol_bin_width);
    fprintf(stderr, "-----\n");
    fprintf(stderr, "min cell density = %.3lf max cell density = %.3lf "
	    "avg cell density = %.3lf units^3\n",
	    stats.min_cell_dense, stats.max_cell_dense, stats.avg_cell_dense);
    fprintf(stderr, "-----\n");
    fprintf(stderr, "cell density histogram:\n");
    fprintf(stderr, "min value\tcount\t\tmax value\n");
    for (k = 0; k < stats.num_dense_bins; k++)
      fprintf(stderr, "%.3lf\t\t%d\t\t%.3lf\n", 
	      stats.min_cell_dense + k * dense_bin_width, stats.dense_hist[k], 
	      stats.min_cell_dense + (k + 1) * dense_bin_width);
    fprintf(stderr, "-------------------------------------------------\n");
  }

}
/*--------------------------------------------------------------------------*/
/*
  aggregates local statistics into global statistics

  nblocks: number of blocks
  vblocks: pointer to array of vblocks
  loc_stats: local statistics
*/
void aggregate_stats(int nblocks, struct vblock_t *vblocks, 
		     struct stats_t *loc_stats) {

  float vol_bin_width; /* width of a volume histogram bin */
  float dense_bin_width; /* width of a density histogram bin */
  struct stats_t glo_stats; /* global stats */
  struct stats_t fin_stats; /* final (global) stats */
  int groupsize; /* MPI usual */
  DIY_Datatype dtype; /* custom datatype */
  MPI_Op op1, op2; /* custom operators */
  int i, j, k;

  MPI_Comm_size(comm, &groupsize);

  if (groupsize > 1) {

    /* create datatype */
    struct map_block_t map[] = {

      { DIY_INT,   OFST, 1, 
	offsetof(struct stats_t, tot_cells)      },
      { DIY_INT,   OFST, 1, 
	offsetof(struct stats_t, tot_faces)      },
      { DIY_INT,   OFST, 1, 
	offsetof(struct stats_t, tot_verts)      },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, avg_cell_verts) },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, avg_cell_faces) },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, avg_face_verts) },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, min_cell_vol)   },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, max_cell_vol)   },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, avg_cell_vol)   },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, min_cell_dense) },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, max_cell_dense) },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, avg_cell_dense) },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, min_cell_time)  },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, max_cell_time)  },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, min_hull_time)  },
      { DIY_FLOAT, OFST, 1, 
	offsetof(struct stats_t, max_hull_time)  },
      { DIY_INT,   OFST, 1, 
	offsetof(struct stats_t, num_vol_bins)   },
      { DIY_INT,   OFST, 1, 
	offsetof(struct stats_t, num_dense_bins) },
      { DIY_INT,   OFST, MAX_HIST_BINS, 
	offsetof(struct stats_t, vol_hist)       },
      { DIY_INT,   OFST, MAX_HIST_BINS, 
	offsetof(struct stats_t, dense_hist)     },

    };

    DIY_Create_struct_datatype(0, 20, map, &dtype);

    MPI_Op_create(&average, 1, &op1);
    MPI_Op_create(&histogram, 1, &op2);

    /* first reduction computes averages and ranges */
    MPI_Reduce(loc_stats, &glo_stats, 1, dtype, op1, 0, comm);
    /* broadcast global stats to all process */
    MPI_Bcast(&glo_stats, 1, dtype, 0, comm);

  }

  else {

    glo_stats.tot_cells      = loc_stats->tot_cells;
    glo_stats.tot_faces      = loc_stats->tot_faces;
    glo_stats.tot_verts      = loc_stats->tot_verts;
    glo_stats.avg_cell_verts = loc_stats->avg_cell_verts;
    glo_stats.avg_cell_faces = loc_stats->avg_cell_faces;
    glo_stats.avg_face_verts = loc_stats->avg_face_verts;
    glo_stats.min_cell_vol   = loc_stats->min_cell_vol;
    glo_stats.max_cell_vol   = loc_stats->max_cell_vol;
    glo_stats.avg_cell_vol   = loc_stats->avg_cell_vol;
    glo_stats.min_cell_dense = loc_stats->min_cell_dense;
    glo_stats.max_cell_dense = loc_stats->max_cell_dense;
    glo_stats.avg_cell_dense = loc_stats->avg_cell_dense;
    glo_stats.min_cell_time  = loc_stats->min_cell_time;
    glo_stats.max_cell_time  = loc_stats->max_cell_time;
    glo_stats.min_hull_time  = loc_stats->min_hull_time;
    glo_stats.max_hull_time  = loc_stats->max_hull_time;
    glo_stats.num_vol_bins   = 50;
    glo_stats.num_dense_bins = 100;

  }

  /* find local cell volume and density histograms */
  vol_bin_width = (glo_stats.max_cell_vol - glo_stats.min_cell_vol) / 
    glo_stats.num_vol_bins; /* volume */
  dense_bin_width = (glo_stats.max_cell_dense - glo_stats.min_cell_dense) / 
    glo_stats.num_dense_bins; /* density */
  for (k = 0; k < glo_stats.num_vol_bins; k++) /* volume */
    glo_stats.vol_hist[k] = 0;
  for (k = 0; k < glo_stats.num_dense_bins; k++) /* density */
    glo_stats.dense_hist[k] = 0;
  for (i = 0; i < nblocks; i++) { /* for all blocks */

    for (j = 0; j < vblocks[i].num_complete_cells; j++) { /* for all cells */

      /* volume */
      for (k = 0; k < glo_stats.num_vol_bins; k++) { /* for all bins */
	if (vblocks[i].vols[j] >= glo_stats.min_cell_vol + k * vol_bin_width && 
	    vblocks[i].vols[j] < 
	    glo_stats.min_cell_vol + (k + 1) * vol_bin_width) {
	  glo_stats.vol_hist[k]++;
	  break;
	}
      } /* for all bins */
      if (k == glo_stats.num_vol_bins)
	glo_stats.vol_hist[k - 1]++; /* catch roundoff error and open
					interval on right side of bin */

      /* density */
      for (k = 0; k < glo_stats.num_dense_bins; k++) { /* for all bins */
	if (vblocks[i].vols[j] > 0.0) {
	  float dense = 1.0 /vblocks[i].vols[j];
	  if (dense >= glo_stats.min_cell_dense + k * dense_bin_width && 
	      dense <  glo_stats.min_cell_dense + (k + 1) * dense_bin_width) {
	    glo_stats.dense_hist[k]++;
	    break;
	  }
	}
      } /* for all bins */
      if (k == glo_stats.num_dense_bins)
	glo_stats.dense_hist[k - 1]++; /* catch roundoff error and open
					interval on right side of bin */

    } /* for all cells */

  } /* for all blocks */

  if (groupsize > 1) {

    /* second reduction computes global histogram */
    MPI_Reduce(&glo_stats, &fin_stats, 1, dtype, op2, 0, comm);

    /* copy global stats back to local stats */
    loc_stats->tot_cells      = fin_stats.tot_cells;
    loc_stats->tot_faces      = fin_stats.tot_faces;
    loc_stats->tot_verts      = fin_stats.tot_verts;
    loc_stats->avg_cell_verts = fin_stats.avg_cell_verts;
    loc_stats->avg_cell_faces = fin_stats.avg_cell_faces;
    loc_stats->avg_face_verts = fin_stats.avg_face_verts;
    loc_stats->min_cell_vol   = fin_stats.min_cell_vol;
    loc_stats->max_cell_vol   = fin_stats.max_cell_vol;
    loc_stats->avg_cell_vol   = fin_stats.avg_cell_vol;
    loc_stats->min_cell_dense = fin_stats.min_cell_dense;
    loc_stats->max_cell_dense = fin_stats.max_cell_dense;
    loc_stats->avg_cell_dense = fin_stats.avg_cell_dense;
    loc_stats->min_cell_time  = fin_stats.min_cell_time;
    loc_stats->max_cell_time  = fin_stats.max_cell_time;
    loc_stats->min_hull_time  = fin_stats.min_hull_time;
    loc_stats->max_hull_time  = fin_stats.max_hull_time;
    loc_stats->num_vol_bins   = fin_stats.num_vol_bins;
    loc_stats->num_dense_bins   = fin_stats.num_dense_bins;
    for (i = 0; i < MAX_HIST_BINS; i++) {
      loc_stats->vol_hist[i] = fin_stats.vol_hist[i];
      loc_stats->dense_hist[i] = fin_stats.dense_hist[i];
    }

    DIY_Destroy_datatype(&dtype);
    MPI_Op_free(&op1);
    MPI_Op_free(&op2);

  }
  else {

    loc_stats->num_vol_bins   = glo_stats.num_vol_bins;
    loc_stats->num_dense_bins = glo_stats.num_dense_bins;
    for (i = 0; i < MAX_HIST_BINS; i++) {
      loc_stats->vol_hist[i] = glo_stats.vol_hist[i];
      loc_stats->dense_hist[i] = glo_stats.dense_hist[i];
    }

  }

}
/*--------------------------------------------------------------------------*/
/*
  reduces averages and ranges

  in: input 1
  inout: input 2 and output
  len: 1
  datatype: unused
*/
void average(void *in, void *inout, int *len, MPI_Datatype *type) {

  /* quiet compiler warnings about unused variables */
  type = type;
  len = len;

  struct stats_t *stats1 = (struct stats_t *)in;
  struct stats_t *stats2 = (struct stats_t *)inout;

  /* weights for weighted averages based on cell counts */
  float w1 = (float)stats1->tot_cells / 
    (stats1->tot_cells + stats2->tot_cells);
  float w2 = (float)stats2->tot_cells / 
    (stats1->tot_cells + stats2->tot_cells);

  /* weighted average of two averages */
  stats2->avg_cell_verts = w1 * stats1->avg_cell_verts +
    w2 * stats2->avg_cell_verts;
  stats2->avg_cell_faces = w1 * stats1->avg_cell_faces +
    w2 * stats2->avg_cell_faces;
  stats2->avg_cell_vol = w1 * stats1->avg_cell_vol +
    w2 * stats2->avg_cell_vol;
  stats2->avg_cell_dense = w1 * stats1->avg_cell_dense +
    w2 * stats2->avg_cell_dense;

  /* new weights for weighted averages based on face counts */
  w1 = (float)stats1->tot_faces / 
    (stats1->tot_faces + stats2->tot_faces);
  w2 = (float)stats2->tot_faces / 
    (stats1->tot_faces + stats2->tot_faces);

  /* weighted average of two averages */
  stats2->avg_face_verts = w1 * stats1->avg_face_verts + 
    w2 * stats2->avg_face_verts;

  stats2->tot_cells += stats1->tot_cells;
  stats2->tot_verts += stats1->tot_verts;

  if (stats1->min_cell_vol < stats2->min_cell_vol)
    stats2->min_cell_vol = stats1->min_cell_vol;
  if (stats1->max_cell_vol > stats2->max_cell_vol)
    stats2->max_cell_vol = stats1->max_cell_vol;
  if (stats1->min_cell_dense < stats2->min_cell_dense)
    stats2->min_cell_dense = stats1->min_cell_dense;
  if (stats1->max_cell_dense > stats2->max_cell_dense)
    stats2->max_cell_dense = stats1->max_cell_dense;
  if (stats1->min_cell_time < stats2->min_cell_time)
    stats2->min_cell_time = stats1->min_cell_time;
  if (stats1->max_cell_time > stats2->max_cell_time)
    stats2->max_cell_time = stats1->max_cell_time;
  if (stats1->min_hull_time < stats2->min_hull_time)
    stats2->min_hull_time = stats1->min_hull_time;
  if (stats1->max_hull_time > stats2->max_hull_time)
    stats2->max_hull_time = stats1->max_hull_time;

  /* ought to do a cross-validation to find correct number of bins
     for now just pick a number */
  stats2->num_vol_bins = 50;
  stats2->num_dense_bins = 100;

}
/*--------------------------------------------------------------------------*/
/*
  reduces histograms

  in: input 1
  inout: input 2 and output
  len: 1
  datatype: unused
*/
void histogram(void *in, void *inout, int *len, MPI_Datatype *type) {

  /* quiet compiler warnings about unused variables */
  type = type;
  len = len;

  struct stats_t *stats1 = (struct stats_t *)in;
  struct stats_t *stats2 = (struct stats_t *)inout;
  int i;

  stats2->tot_cells      = stats1->tot_cells;
  stats2->tot_faces      = stats1->tot_faces;
  stats2->tot_verts      = stats1->tot_verts;
  stats2->avg_cell_verts = stats1->avg_cell_verts;
  stats2->avg_cell_faces = stats1->avg_cell_faces;
  stats2->avg_face_verts = stats1->avg_face_verts;
  stats2->min_cell_vol   = stats1->min_cell_vol;
  stats2->max_cell_vol   = stats1->max_cell_vol;
  stats2->avg_cell_vol   = stats1->avg_cell_vol;
  stats2->min_cell_dense = stats1->min_cell_dense;
  stats2->max_cell_dense = stats1->max_cell_dense;
  stats2->avg_cell_dense = stats1->avg_cell_dense;
  stats2->min_cell_time  = stats1->min_cell_time;
  stats2->max_cell_time  = stats1->max_cell_time;
  stats2->min_hull_time  = stats1->min_hull_time;
  stats2->max_hull_time  = stats1->max_hull_time;
  stats2->num_vol_bins   = stats1->num_vol_bins;
  stats2->num_dense_bins = stats1->num_dense_bins;
  for (i = 0; i < stats2->num_vol_bins; i++)
    stats2->vol_hist[i] += stats1->vol_hist[i];
  for (i = 0; i < stats2->num_dense_bins; i++)
    stats2->dense_hist[i] += stats1->dense_hist[i];

}
/*--------------------------------------------------------------------------*/
/*
  prepare for output

  nblocks: number of blocks
  vblocks: pointer to array of vblocks
  hdrs: block headers
*/
void prep_out(int nblocks, struct vblock_t *vblocks, int **hdrs) {

  struct bb_t bounds; /* block bounds */
  int i, j;

  /* save extents */
  for (i = 0; i < nblocks; i++) {
    DIY_Block_bounds(0, i, &bounds);
    vblocks[i].mins[0] = bounds.min[0];
    vblocks[i].mins[1] = bounds.min[1];
    vblocks[i].mins[2] = bounds.min[2];
    vblocks[i].maxs[0] = bounds.max[0];
    vblocks[i].maxs[1] = bounds.max[1];
    vblocks[i].maxs[2] = bounds.max[2];
  }

  /* save vertices (float version) */
  for (i = 0; i < nblocks; i++) {
    vblocks[i].save_verts = (float *)malloc(vblocks[i].num_verts * 3 * 
					       sizeof(float));
    for (j = 0; j < vblocks[i].num_verts; j++) {
      /* added by jingyuan
       * check and remove vertices that are out of data extent
       */
      if (vblocks[i].verts[3 * j    ] < data_mins[0] ||
          vblocks[i].verts[3 * j    ] > data_maxs[0] ||
          vblocks[i].verts[3 * j + 1] < data_mins[1] ||
          vblocks[i].verts[3 * j + 1] > data_maxs[1] ||
          vblocks[i].verts[3 * j + 2] < data_mins[2] ||
          vblocks[i].verts[3 * j + 2] > data_maxs[2]) {
        vblocks[i].save_verts[3 * j    ] = 0;
        vblocks[i].save_verts[3 * j + 1] = 0;
        vblocks[i].save_verts[3 * j + 2] = 0;
      }
      else {
        vblocks[i].save_verts[3 * j]     = vblocks[i].verts[3 * j];
        vblocks[i].save_verts[3 * j + 1] = vblocks[i].verts[3 * j + 1];
        vblocks[i].save_verts[3 * j + 2] = vblocks[i].verts[3 * j + 2];
      }
      /* end of changes */
    }
    free (vblocks[i].verts);
    vblocks[i].verts = NULL;
  }

  /* save headers */
  for (i = 0; i < nblocks; i++) {

    hdrs[i][NUM_VERTS] = vblocks[i].num_verts;
    hdrs[i][NUM_CELLS] = vblocks[i].num_cells;
    hdrs[i][TOT_NUM_CELL_VERTS] = vblocks[i].tot_num_cell_verts;
    hdrs[i][NUM_COMPLETE_CELLS] = vblocks[i].num_complete_cells;
    hdrs[i][TOT_NUM_CELL_FACES] = vblocks[i].tot_num_cell_faces;
    hdrs[i][TOT_NUM_FACE_VERTS] = vblocks[i].tot_num_face_verts;
    hdrs[i][NUM_ORIG_PARTICLES] = vblocks[i].num_orig_particles;

  }

}
/*--------------------------------------------------------------------------*/
/*
  creates and initializes blocks and headers

  num_blocks: number of blocks
  vblocks: pointer to array of vblocks
  hdrs: pointer to array of headers

  side effects: allocates memory for blocks and headers
*/
void create_blocks(int num_blocks, struct vblock_t **vblocks, int ***hdrs) {

  int i, j;

  /* allocate blocks and headers */
  *vblocks = (struct vblock_t*)malloc(sizeof(struct vblock_t) * 
				      num_blocks);
  *hdrs = (int **)malloc(sizeof(int*) * num_blocks);

  for (i = 0; i < num_blocks; i++) {

    (*vblocks)[i].num_verts = 0;
    (*vblocks)[i].num_cells = 0;
    (*vblocks)[i].verts = NULL;
    (*vblocks)[i].save_verts = NULL;
    (*vblocks)[i].num_cell_verts = NULL;
    (*vblocks)[i].tot_num_cell_verts = 0;
    (*vblocks)[i].cells = NULL;
    (*vblocks)[i].sites = NULL;
    (*vblocks)[i].temp_num_complete_cells = 0;
    (*vblocks)[i].temp_complete_cells = NULL;
    (*vblocks)[i].num_complete_cells = 0;
    (*vblocks)[i].complete_cells = NULL;
    (*vblocks)[i].areas = NULL;
    (*vblocks)[i].vols = NULL;
    (*vblocks)[i].tot_num_cell_faces = 0;
    (*vblocks)[i].tot_num_face_verts = 0;
    (*vblocks)[i].num_cell_faces = NULL;
    (*vblocks)[i].num_face_verts = NULL;
    (*vblocks)[i].face_verts = NULL;

    (*hdrs)[i] = (int *)malloc(sizeof(int) * DIY_MAX_HDR_ELEMENTS);
    for (j = 0; j < DIY_MAX_HDR_ELEMENTS; j++)
      ((*hdrs)[i])[j] = 0;

  }

}
/*---------------------------------------------------------------------------*/
/*
  frees blocks and headers

  num_blocks: number of blocks
  vblocks: pointer to array of vblocks
  hdrs: pointer to array of headers
*/
void destroy_blocks(int num_blocks, struct vblock_t *vblocks, int **hdrs) {

  int i;

  for (i = 0; i < num_blocks; i++) {
    if (hdrs[i])
      free(hdrs[i]);
    if (vblocks[i].verts)
      free(vblocks[i].verts);
    if (vblocks[i].save_verts)
      free(vblocks[i].save_verts);
    if (vblocks[i].num_cell_verts)
      free(vblocks[i].num_cell_verts);
    if (vblocks[i].cells)
      free(vblocks[i].cells);
    if (vblocks[i].sites)
      free(vblocks[i].sites);
    if (vblocks[i].temp_complete_cells)
      free(vblocks[i].temp_complete_cells);
    if (vblocks[i].complete_cells)
      free(vblocks[i].complete_cells);
    if (vblocks[i].areas)
      free(vblocks[i].areas);
    if (vblocks[i].vols)
      free(vblocks[i].vols);
    if (vblocks[i].num_cell_faces)
      free(vblocks[i].num_cell_faces);
    if (vblocks[i].num_face_verts)
      free(vblocks[i].num_face_verts);
    if (vblocks[i].face_verts)
      free(vblocks[i].face_verts);
  }
  free(hdrs);
  free(vblocks);

}
/*---------------------------------------------------------------------------*/
/*
  determnines complete cells: cells that don't contain qhull'sinfinite vertex or
  any other vertices outside of the block bounds

  vblock: one voronoi block
  lid: local id of block
  ghost: extra size per side beyond block bounds to consider vertices as valid

  side effects: allocates memory for complete cells in voronoi block
*/
void complete_cells(struct vblock_t *vblock, int lid, float ghost) {

  struct bb_t bounds; /* block bounds */
  int vid, vid1; /* vertex id */
  int j, k, n, m;
  double d2_min = 0.0; /* dia^2 of circumscribing sphere of volume min_vol */
  int start_n; /* index into cells at start of new cell */
  int too_small; /* whether cell volume is definitely below threshold */

  /* allocate memory based on number of cells and average number of
     faces and vertices in a cell
     this is wasteful if filtering on volume because we will probably
     only need a fraction of this memory
     todo: fix this with my own memory manager */
  vblock->temp_complete_cells =
	  (int *)malloc(vblock->num_cells * sizeof(int));
  vblock->complete_cells =
	  (int *)malloc(vblock->num_cells * sizeof(int));
  vblock->areas =
    (float *)malloc(vblock->num_cells * sizeof(float));
  vblock->vols =
    (float *)malloc(vblock->num_cells * sizeof(float));
  vblock->num_cell_faces =
    (int *)malloc(vblock->num_cells * sizeof(int));
  vblock->num_face_verts =
	  (int *)malloc(vblock->num_cells * AVG_CELL_FACES * sizeof(int));
  vblock->face_verts =
	  (int *)malloc(vblock->num_cells * AVG_CELL_FACES *
			AVG_FACE_VERTS * sizeof(int));

  /* init */
  memset(vblock->num_face_verts, 0,
	 vblock->num_cells * AVG_CELL_FACES * sizeof(int));
  memset(vblock->face_verts, 0, vblock->num_cells * AVG_CELL_FACES *
	 AVG_FACE_VERTS * sizeof(int));
  vblock->tot_num_cell_faces = 0;
  vblock->tot_num_face_verts = 0;

  DIY_Block_bounds(0, lid, &bounds);

  /* minimum cell diameter squared */
  if (min_vol > 0.0)
    /* d^2 = 4 * (3/4 * min_vol / pi)^ (2/3) */
    d2_min = 1.539339 * pow(min_vol, 0.66667);

  /* find complete cells */

  vblock->temp_num_complete_cells = 0;
  n = 0; /* index into vblock->cells */

  /* for all cells */
  for (j = 0; j < vblock->num_cells; j++) {

    /* init */
    if (!vblock->num_cell_verts[j])
      continue;

    int complete = 1;
    too_small = (min_vol > 0.0 ? 1 : 0);

    /* for all vertex indices in the current cell */
    for (k = 0; k < vblock->num_cell_verts[j]; k++) {

      if (k == 0)
	start_n = n;

      vid = vblock->cells[n];

      if ( /* vertex can fail for the following reasons */
	  /* qhull's "infinite vertex" */
	  vblock->cells[n] == 0 ||
	  /* out of block bounds */
	  vblock->verts[3 * vid]     < bounds.min[0] - ghost ||
	  vblock->verts[3 * vid]     > bounds.max[0] + ghost ||
	  vblock->verts[3 * vid + 1] < bounds.min[1] - ghost ||
	  vblock->verts[3 * vid + 1] > bounds.max[1] + ghost ||
	  vblock->verts[3 * vid + 2] < bounds.min[2] - ghost ||
	  vblock->verts[3 * vid + 2] > bounds.max[2] + ghost ||
	  /* out of overall data bounds */
	  vblock->verts[3 * vid]     < data_mins[0] ||
	  vblock->verts[3 * vid]     > data_maxs[0] ||
	  vblock->verts[3 * vid + 1] < data_mins[1] ||
	  vblock->verts[3 * vid + 1] > data_maxs[1] ||
	  vblock->verts[3 * vid + 2] < data_mins[2] ||
	  vblock->verts[3 * vid + 2] > data_maxs[2]
	   ) {
	complete = 0;
	n += (vblock->num_cell_verts[j] - k); /* skip rest of this cell */
	break;
      } /* if */

      /* check minimum volume if enabled and it has not been excceded yet */
      if (too_small) {

	/* for all vertices in this cell */
	for (m = start_n; m < start_n + vblock->num_cell_verts[j]; m++) {
	  vid1 = vblock->cells[m];
	  double d2 =
	    (vblock->verts[3 * vid] - vblock->verts[3 * vid1]) *
	    (vblock->verts[3 * vid] - vblock->verts[3 * vid1]) +
	    (vblock->verts[3 * vid + 1] - vblock->verts[3 * vid1 + 1]) *
	    (vblock->verts[3 * vid + 1] - vblock->verts[3 * vid1 + 1]) +
	    (vblock->verts[3 * vid + 2] - vblock->verts[3 * vid1 + 2]) *
	    (vblock->verts[3 * vid + 2] - vblock->verts[3 * vid1 + 2]);
	  if (d2 > d2_min) {
	    too_small = 0;
	    break;
	  }

	} /* all vertices in this cell */

      } /* samll volume threshold */

      /* check if volume is too small at the end of the cell */
      if (k == vblock->num_cell_verts[j] - 1 && too_small)
	complete = 0;

      n++;

    } /* for all vertex indices in this cell */

    if (complete)
      (vblock->temp_complete_cells)[vblock->temp_num_complete_cells++] = j;

  } /* for all cells */

}
/*--------------------------------------------------------------------------*/
/*
  prepares vertices for convex hull computation

  vblock: one voronoi block
  cell: the current complete cell
  vertices: vertex array, allocated by caller
  vmap: vertex map from these vertices back to original vertex ids, allocated
  by caller
*/
void prep_vertices(struct vblock_t *vblock, int cell, double *vertices,
		   int *vmap) {

  int i, k;
  int num_verts = vblock->num_cell_verts[vblock->temp_complete_cells[cell]];

  /* compute index of starting vertex */
  int start_vert = 0; /* index of first vertex */
  for (i = 0; i < (vblock->temp_complete_cells)[cell]; i++)
    start_vert += (vblock->num_cell_verts)[i];

  /* copy vertices */
  for (k = 0; k < num_verts; k++) {
    int vid = (vblock->cells)[start_vert + k]; /* current vertex id */
    vertices[3 * k]     = (vblock->verts)[3 * vid];
    vertices[3 * k + 1] = (vblock->verts)[3 * vid + 1];
    vertices[3 * k + 2] = (vblock->verts)[3 * vid + 2];
    vmap[k] = vid;

  }

}
/*--------------------------------------------------------------------------*/
/*
  generates test particles for a  block

  lid: local id of block
  particles: pointer to particle vector in this order: 
  particle0x, particle0y, particle0z, particle1x, particle1y, particle1z, ...
  jitter: maximum amount to randomly move particles

  returns: number of particles in this block

  side effects: allocates memory for particles, caller's responsibility to free
*/
int gen_particles(int lid, float **particles, float jitter) {

  int sizes[3]; /* number of grid points */
  int i, j, k;
  int n = 0;
  int num_particles;
  float jit; /* random jitter amount, 0 - MAX_JITTER */

  /* linear positions */

  /* allocate particles */
  struct bb_t bounds;
  DIY_Block_bounds(0, lid, &bounds);
  sizes[0] = (int)(bounds.max[0] - bounds.min[0] + 1);
  sizes[1] = (int)(bounds.max[1] - bounds.min[1] + 1);
  sizes[2] = (int)(bounds.max[2] - bounds.min[2] + 1);
  num_particles = sizes[0] * sizes[1] * sizes[2];
  *particles = (float *)malloc(num_particles * 3 * sizeof(float));

  /* assign particles */

  n = 0;
  for (i = 0; i < sizes[0]; i++) {
    for (j = 0; j < sizes[1]; j++) {
      for (k = 0; k < sizes[2]; k++) {

	// start with particles on a grid
	(*particles)[3 * n] = bounds.min[0] + i;
	(*particles)[3 * n + 1] = bounds.min[1] + j;
	(*particles)[3 * n + 2] = bounds.min[2] + k;

	// and now jitter them
	jit = rand() / (float)RAND_MAX * 2 * jitter - jitter;
	if ((*particles)[3 * n] - jit >= bounds.min[0] &&
	    (*particles)[3 * n] - jit <= bounds.max[0])
	  (*particles)[3 * n] -= jit;
	else if ((*particles)[3 * n] + jit >= bounds.min[0] &&
		 (*particles)[3 * n] + jit <= bounds.max[0])
	  (*particles)[3 * n] += jit;

	jit = rand() / (float)RAND_MAX * 2 * jitter - jitter;
	if ((*particles)[3 * n + 1] - jit >= bounds.min[1] &&
	    (*particles)[3 * n + 1] - jit <= bounds.max[1])
	  (*particles)[3 * n + 1] -= jit;
	else if ((*particles)[3 * n + 1] + jit >= bounds.min[1] &&
		 (*particles)[3 * n + 1] + jit <= bounds.max[1])
	  (*particles)[3 * n + 1] += jit;

	jit = rand() / (float)RAND_MAX * 2 * jitter - jitter;
	if ((*particles)[3 * n + 2] - jit >= bounds.min[2] &&
	    (*particles)[3 * n + 2] - jit <= bounds.max[2])
	  (*particles)[3 * n + 2] -= jit;
	else if ((*particles)[3 * n + 2] + jit >= bounds.min[2] &&
		 (*particles)[3 * n + 2] + jit <= bounds.max[2])
	  (*particles)[3 * n + 2] += jit;


	n++;

      }
    }
  }

  /* quadratic positions */

/*   /\* allocate particles *\/ */
/*   struct bb_t bounds; */
/*   DIY_Block_bounds(0, lid, &bounds); */
/*   sizes[0] = (int)sqrtf(bounds.max[0] - bounds.min[0]) + 1; */
/*   sizes[1] = (int)sqrtf(bounds.max[1] - bounds.min[1]) + 1; */
/*   sizes[2] = (int)sqrtf(bounds.max[2] - bounds.min[2]) + 1; */
/*   num_particles = sizes[0] * sizes[1] * sizes[2]; */
/*   *particles = (float *)malloc(num_particles * 3 * sizeof(float)); */

/*   /\* assign particles *\/ */
/*   n = 0; */
/*   for (i = 0; i < sizes[0]; i++) { */
/*     for (j = 0; j < sizes[1]; j++) { */
/*       for (k = 0; k < sizes[2]; k++) { */

/* 	// start with particles on a grid */
/* 	(*particles)[3 * n] = bounds.min[0] + i * i; */
/* 	(*particles)[3 * n + 1] = bounds.min[1] + j * j; */
/* 	(*particles)[3 * n + 2] = bounds.min[2] + k * k; */

/* 	// and now jitter them */
/* 	jit = rand() / (float)RAND_MAX * 2 * jitter - jitter; */
/* 	if ((*particles)[3 * n] - jit >= bounds.min[0] && */
/* 	    (*particles)[3 * n] - jit <= bounds.max[0]) */
/* 	  (*particles)[3 * n] -= jit; */
/* 	else if ((*particles)[3 * n] + jit >= bounds.min[0] && */
/* 		 (*particles)[3 * n] + jit <= bounds.max[0]) */
/* 	  (*particles)[3 * n] += jit; */

/* 	jit = rand() / (float)RAND_MAX * 2 * jitter - jitter; */
/* 	if ((*particles)[3 * n + 1] - jit >= bounds.min[1] && */
/* 	    (*particles)[3 * n + 1] - jit <= bounds.max[1]) */
/* 	  (*particles)[3 * n + 1] -= jit; */
/* 	else if ((*particles)[3 * n + 1] + jit >= bounds.min[1] && */
/* 		 (*particles)[3 * n + 1] + jit <= bounds.max[1]) */
/* 	  (*particles)[3 * n + 1] += jit; */

/* 	jit = rand() / (float)RAND_MAX * 2 * jitter - jitter; */
/* 	if ((*particles)[3 * n + 2] - jit >= bounds.min[2] && */
/* 	    (*particles)[3 * n + 2] - jit <= bounds.max[2]) */
/* 	  (*particles)[3 * n + 2] -= jit; */
/* 	else if ((*particles)[3 * n + 2] + jit >= bounds.min[2] && */
/* 		 (*particles)[3 * n + 2] + jit <= bounds.max[2]) */
/* 	  (*particles)[3 * n + 2] += jit; */


/* 	n++; */

/*       } */
/*     } */
/*   } */

  /* debug: print particles */
/*   for (n = 0; n < num_particles; n++) */
/*     fprintf(stderr, "block = %d particle[%d] = [%.2lf %.2lf %.2lf]\n", */
/* 	    lid, n, (*particles)[3 * n], (*particles)[3 * n + 1], */
/* 	    (*particles)[3 * n + 2]); */

  return num_particles;

}
/*--------------------------------------------------------------------------*/
/*
  generates voronoi output from qhull

  facetlist: qhull list of convex hull facets
  vblock: pointer to one voronoi block, allocated by caller
  num_particles: number of particles used to generate the tessellation
  side effects: allocates data structures inside of vblock, caller's
  responsibility to free
*/
void gen_voronoi_output(facetT *facetlist, struct vblock_t *vblock, 
			int num_particles) {

  int i, numcenters, numvertices= 0, numneighbors, numinf, vertex_i, vertex_n;
  facetT *facet, *neighbor, **neighborp;
  setT *vertices;
  vertexT *vertex;
  boolT isLower;
  unsigned int numfacets= (unsigned int) qh num_facets;

  /* init, get counts */
  vertices = qh_markvoronoi(facetlist, NULL, 0, &isLower, &numcenters);
  FOREACHvertex_i_(vertices) {
    if (vertex) {
      numvertices++;
      numneighbors = numinf = 0;
      FOREACHneighbor_(vertex) {
        if (neighbor->visitid == 0)
          numinf= 1;
        else if (neighbor->visitid < numfacets)
          numneighbors++;
      }
      if (numinf && !numneighbors) {
        SETelem_(vertices, vertex_i)= NULL;
        numvertices--;
      }
    }
  }
  /* number of verts and cells may appear to be reversed, but this is
     qhull's nomenclature and is actually correct
     also, num_cells is temporary because it includes cells formed by 
     ghost perticles
     eventually vblock->num_cells will be (permanentlly) set to only the 
     cells derived from original particles
  */
  vblock->num_verts = numcenters;
  int temp_num_cells = numvertices;

  /* vertices */
  vblock->verts = (double *)malloc(sizeof(double) * 3 * vblock->num_verts);
  vblock->verts[0] = vblock->verts[1] = vblock->verts[2] = qh_INFINITE;
  i = 1; /* already did the infinity vertex, index 0 */
  FORALLfacet_(facetlist) {
    if (facet->visitid && facet->visitid < numfacets) {
      if (!facet->center)
	facet->center = qh_facetcenter(facet->vertices);
      vblock->verts[3 * i] = facet->center[0];
      vblock->verts[3 * i + 1] = facet->center[1];
      vblock->verts[3 * i + 2] = facet->center[2];
      i++;
    }
  }

  /* number of vertices in each cell. temp_num_cells <= num_particles, 
     allocate to quantity num_particles to be safe */
  vblock->num_cell_verts = (int *)malloc(sizeof(int) * num_particles);
  memset(vblock->num_cell_verts, 0, sizeof(int) * num_particles);
  i = 0;
  FOREACHvertex_i_(vertices) {
    numneighbors = 0;
    numinf = 0;
    if (vertex) {
      FOREACHneighbor_(vertex) {
        if (neighbor->visitid == 0)
          numinf= 1;
        else if (neighbor->visitid < numfacets)
          numneighbors++;
      }
      if (numinf)
	numneighbors++;
      vblock->num_cell_verts[i++] = numneighbors;
    }
  }

  /* allocate the cell vertices */
  vblock->tot_num_cell_verts = 0;
  for (i = 0; i < temp_num_cells; i++)
    vblock->tot_num_cell_verts += vblock->num_cell_verts[i];
  vblock->cells = (int *)malloc(sizeof(int) * vblock->tot_num_cell_verts);

  /* cell vertices */
  i = 0;
  FOREACHvertex_i_(vertices) {
    if (vertex) {
      numinf = 0;
      FOREACHneighbor_(vertex) {
	if (neighbor->visitid < numfacets) {
	  if (!numinf || neighbor->visitid > 0) {
	    vblock->cells[i++] = neighbor->visitid;
	    if (neighbor->visitid == 0)
	      numinf++;
	  }
	}
	else if (numinf && neighbor->visitid < numfacets)
	  vblock->cells[i++] = neighbor->visitid;
      }
    }
  }


  /* debug: print output */
/*   fprintf(stderr, "%d cells:\n", vblock->num_cells);  */
/*   int n = 0; */
/*   for(i = 0; i < vblock->num_cells; i++) { */
/*   fprintf(stderr, "site %.1lf %.1lf %.1lf: %d vertices: ", vblock->sites[3 * i], */
/* 	  vblock->sites[3 * i + 1], vblock->sites[3 * i + 2],  */
/* 	  vblock->num_cell_verts[i]); */
/*     for (j = 0; j < vblock->num_cell_verts[i]; j++) */
/*       fprintf(stderr, "%d ", vblock->cells[n++]); */
/*     fprintf(stderr, "\n"); */
/*   } */

  /* clean up */
  qh_settempfree(&vertices);

}
/*--------------------------------------------------------------------------*/
/*
  generates convex hull output

  facetlist: qhull list of convex hull facets
  cblock: pointer to one convex hull block

  returns: 0 = hull did not pass volume filter test, discard
  1 = hull passed volume filter test, save

  side effects: if return value is 1, allocates num_cell_faces, 
  num_face_verts, face_verts, areas,
  vols data structures inside of vblock, caller's responsibility to free
*/
int gen_convex_output(facetT *facetlist, struct cblock_t *cblock) {

  facetT *facet;
  setT *vertices;
  vertexT *vertex, **vertexp;
  int num_verts; /* number of vertices in current face */
  int face; /* current face */
  int vert; /* current vertex */

  if ((min_vol > 0 && qh totvol < min_vol) ||
      (max_vol > 0 && qh totvol > max_vol))
    return 0;

  /* get number of faces */
  cblock->num_cell_faces = 0;
  FORALLfacet_(facetlist)
    cblock->num_cell_faces++;

  /* allocate memory and assign values */
  cblock->num_face_verts = 
	  (int *)malloc(cblock->num_cell_faces * sizeof(int));
  cblock->face_verts = 
	  (int **)malloc(cblock->num_cell_faces * sizeof(int *));

  face = 0;
  FORALLfacet_(facetlist) {
    vertices= qh_facet3vertex(facet);
    num_verts = qh_setsize(vertices);
    cblock->num_face_verts[face] = num_verts;
    cblock->face_verts[face] = (int *)malloc(num_verts * sizeof(int));
    vert = 0;
    FOREACHvertex_(vertices) {
      cblock->face_verts[face][vert] = qh_pointid(vertex->point);
      vert++;
    }
    qh_settempfree(&vertices);
    face++;
  }

  /* area and volume */
  cblock->area = qh totarea;
  cblock->vol = qh totvol;

  return 1;

}
/*--------------------------------------------------------------------------*/
/*
  copies convex hull output to voronoi block

  cblock: pointer to one convex hull block
  vblock: pointer to one voronoi block, allocated by caller
  vmap: map of convex hull vertex ids back to original voronoi ids
  cell: index of current cell

  side effects: allocates complete_cells, areas, vols, num_cell_faces
  data structures inside of vblock, caller's responsibility to free
*/
void convex_to_voronoi(struct cblock_t *cblock, struct vblock_t *vblock,
		       int *vmap, int cell) {

  int j, k;

  /* todo: replace the following assertions with dynamic reallocation
     of vblock->num_face_verts and vblock->face_verts */

  /* copy */
  (vblock->areas)[cell] = cblock->area;
  (vblock->vols)[cell] = cblock->vol;
  (vblock->num_cell_faces)[cell] = cblock->num_cell_faces;

  /* confirm that the size of vblock->num_face_verts is large enough */
  assert(vblock->tot_num_cell_faces + (vblock->num_cell_faces)[cell] <=
	 vblock->temp_num_complete_cells * AVG_CELL_FACES);

  for (j = 0; j < (vblock->num_cell_faces)[cell]; j++) {

    (vblock->num_face_verts)[vblock->tot_num_cell_faces] = 
      (cblock->num_face_verts)[j];

    /* confirm that the size of vblock->face_verts is large enough */
    assert(vblock->tot_num_face_verts + 
	   (vblock->num_face_verts)[vblock->tot_num_cell_faces] <=
	   vblock->temp_num_complete_cells * AVG_CELL_FACES * AVG_FACE_VERTS);

    for (k = 0; k < (vblock->num_face_verts)[vblock->tot_num_cell_faces]; k++) {
      (vblock->face_verts)[vblock->tot_num_face_verts] = 
	vmap[(cblock->face_verts)[j][k]];
      vblock->tot_num_face_verts++;

    }

    vblock->tot_num_cell_faces++;

  }

  /* free convex hull memory */
  for (j = 0; j < cblock->num_cell_faces; j++)
    free((cblock->face_verts)[j]);
  free(cblock->face_verts);
  free(cblock->num_face_verts);

}
/*--------------------------------------------------------------------------*/
/*
  prints a block

  vblock: current voronoi block
  gid: global block id
*/
void print_block(struct vblock_t *vblock, int gid) {

  int i;

  fprintf(stderr, "block gid = %d, %d complete cells: ", 
	  gid, vblock->num_complete_cells);
  for (i = 0; i < vblock->num_complete_cells; i++)
    fprintf(stderr, "%d ", vblock->complete_cells[i]);
  fprintf(stderr, "\n");

}
/*--------------------------------------------------------------------------*/
/*
  prints particles

  prticles: particle array
  num_particles: number of particles
  gid: block global id
*/
void print_particles(float *particles, int num_particles, int gid) {

  int n;

  for (n = 0; n < num_particles; n++)
    fprintf(stderr, "block = %d particle[%d] = [%.1lf %.1lf %.1lf]\n",
	    gid, n, particles[3 * n], particles[3 * n + 1],
	    particles[3 * n + 2]);

}
/*--------------------------------------------------------------------------*/
/*
  transforms particles for enqueueing to wraparound neighbors
  p: pointer to particle
  wrap_dir: wrapping direcion
*/
void transform_particle(char *p, unsigned char wrap_dir) {

  /* debug */
  float particle[3]; /* original particle */
  particle[0] = ((float *)p)[0];
  particle[1] = ((float *)p)[1];
  particle[2] = ((float *)p)[2];

  int print = 0; /* whether to print the output */

  /* wrapping toward the left transforms to the right */
  if ((wrap_dir & DIY_X0) == DIY_X0)  {
    ((float *)p)[0] += (data_maxs[0] - data_mins[0]);
    print = 1;
  }
  /* and vice versa */
  if ((wrap_dir & DIY_X1) == DIY_X1)  {
    ((float *)p)[0] -= (data_maxs[0] - data_mins[0]);
    print = 1;
  }

  /* similar for y, z */
  if ((wrap_dir & DIY_Y0) == DIY_Y0)  {
    ((float *)p)[1] += (data_maxs[1] - data_mins[1]);
    print = 1;
  }
  if ((wrap_dir & DIY_Y1) == DIY_Y1)  {
    ((float *)p)[1] -= (data_maxs[1] - data_mins[1]);
    print = 1;
  }

  if ((wrap_dir & DIY_Z0) == DIY_Z0)  {
    ((float *)p)[2] += (data_maxs[2] - data_mins[2]);
    print = 1;
  }
  if ((wrap_dir & DIY_Z1) == DIY_Z1)  {
    ((float *)p)[2] -= (data_maxs[2] - data_mins[2]);
    print = 1;
  }

  /* debug */
/*   if (print) */
/*     fprintf(stderr, "before miror [%.2lf %.2lf %.2lf] " */
/* 	    "after miror [%.2lf %.2lf %.2lf]\n", */
/* 	    particle[0], particle[1], particle[2], */
/* 	    ((float *)p)[0], ((float *)p)[1], ((float *)p)[2]); */

/*   else */
/*     fprintf(stderr, "no miror [%.2lf %.2lf %.2lf]\n", */
/* 	    ((float *)p)[0], ((float *)p)[1], ((float *)p)[2]); */

}
/*--------------------------------------------------------------------------*/
