/*---------------------------------------------------------------------------
 *
 * parallel netcdf I/O for voronoi tesselation
 *
 * Wei-keng Liao (Northwestern University)
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

#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <assert.h>
#include <mpi.h>

#ifdef USEPNETCDF
#include <pnetcdf.h>
#endif

#include "voronoi.h"
#include "diy.h"
#include "io.h"

/*---------------------------------------------------------------------------
 *  pnetcdf file schema
 *
 *      dimensions:
 *              num_g_blocks; ie, tot_blocks
 *              XYZ = 3;
 *              num_g_verts;
 *              num_g_complete_cells;
 *              tot_num_g_cell_faces;
 *              tot_num_g_face_verts;
 *              num_g_orig_particles;
 *              num_g_neighbors;
 *      variables:
 *              int num_verts(num_g_blocks) ;
 *              int num_complete_cells(num_g_blocks) ;
 *              int tot_num_cell_faces(num_g_blocks) ;
 *              int tot_num_face_verts(num_g_blocks) ;
 *              int num_orig_particles(num_g_blocks) ;
 *              int64 block_off_num_verts(num_g_blocks) ;
 *              int64 block_off_num_complete_cells(num_g_blocks) ;
 *              int64 block_off_tot_num_cell_faces(num_g_blocks) ;
 *              int64 block_off_tot_num_face_verts(num_g_blocks) ;
 *              int64 block_off_num_orig_particles(num_g_blocks) ;
 *              float mins(tot_blocks, XYZ) ;
 *              float maxs(tot_blocks, XYZ) ;
 *              float save_verts(num_g_verts, XYZ) ;
 *              float sites(num_g_orig_particles) ;
 *              int complete_cells(num_g_complete_cells) ;
 *              float areas(num_g_complete_cells) ;
 *              float vols(num_g_complete_cells) ;
 *              int num_cell_faces(num_g_complete_cells) ;
 *              int num_face_verts(tot_num_g_cell_faces) ;
 *              int face_verts(tot_num_g_face_verts) ;
 *              int neighbors(num_g_neighbors) ;
 *              int g_block_ids(tot_blocks) ;
 *
---------------------------------------------------------------------------*/
/*
  writes output in diy format

  nblocks: number of blocks
  vblocks: pointer to array of vblocks
  hdrs: block headers
  out_file: output file name
*/
void diy_write(int nblocks, struct vblock_t *vblocks, int ** hdrs,
         char *out_file) {

  int i;

  /* pointers to voronoi blocks, needed for writing output */
  void **pvblocks;
  pvblocks = (void**)malloc(sizeof(void*) * nblocks);
  for (i = 0; i < nblocks; i++)
    pvblocks[i] = &vblocks[i];

  /* write output */
  DIY_Write_open_all(0, out_file, 0); /* uncompressed for now */
  DIY_Write_blocks_all(0, pvblocks, nblocks, hdrs, &create_datatype);
  DIY_Write_close_all(0);

  free(pvblocks);

}
/*--------------------------------------------------------------------------*/
/*
  writes output in pnetcdf format

  nblocks: local number of blocks
  vblocks: pointer to array of vblocks
  out_file: output file name
  comm: MPI communicator
*/
void pnetcdf_write(int nblocks, struct vblock_t *vblocks,
       char *out_file, MPI_Comm comm) {

#ifdef USEPNETCDF
  int err;
  int ncid, cmode, varids[23], dimids[8], dimids_2D[2];
  MPI_Offset start[2], count[2];

  MPI_Offset quants[NUM_QUANTS]; /* quantities per block */
  MPI_Offset proc_quants[NUM_QUANTS]; /* quantities per process */
  MPI_Offset tot_quants[NUM_QUANTS]; /* total quantities all global blocks */
  MPI_Offset block_ofsts[NUM_QUANTS]; /* starting offsets for each block */

  /* init */
  int i;
  for (i = 0; i < NUM_QUANTS; i++) {
    quants[i] = 0;
    proc_quants[i] = 0;
    tot_quants[i] = 0;
    block_ofsts[i] = 0;
  }

  /* sum quantities over local blocks */
  int b;
  for (b = 0; b < nblocks; b++) {
    proc_quants[NUM_VERTS] += vblocks[b].num_verts;
    proc_quants[NUM_COMP_CELLS] += vblocks[b].num_complete_cells;
    proc_quants[NUM_CELL_FACES] += vblocks[b].tot_num_cell_faces;
    proc_quants[NUM_FACE_VERTS] += vblocks[b].tot_num_face_verts;
    proc_quants[NUM_ORIG_PARTS] += vblocks[b].num_orig_particles;
    proc_quants[NUM_NEIGHBORS] += DIY_Num_neighbors(0, b);
  }
  proc_quants[NUM_BLOCKS] = nblocks;

  /* sum per process values to be global ones */
  MPI_Allreduce(proc_quants, tot_quants, NUM_QUANTS, MPI_OFFSET, MPI_SUM, comm);

  /* prefix sum proc offsets */
  MPI_Exscan(proc_quants, &block_ofsts, NUM_QUANTS, MPI_OFFSET, MPI_SUM, comm);

  /* create a new file for writing */
  cmode = NC_CLOBBER | NC_64BIT_DATA;
  err = ncmpi_create(comm, out_file, cmode, MPI_INFO_NULL, &ncid); ERR;

  /* define dimensions */
  err = ncmpi_def_dim(ncid, "num_g_blocks", tot_quants[NUM_BLOCKS],
          &dimids[0]); ERR;
  err = ncmpi_def_dim(ncid, "XYZ", 3, &dimids[1]); ERR;
  err = ncmpi_def_dim(ncid, "num_g_verts", tot_quants[NUM_VERTS],
          &dimids[2]); ERR;
  err = ncmpi_def_dim(ncid, "num_g_complete_cells", tot_quants[NUM_COMP_CELLS],
          &dimids[3]); ERR;
  err = ncmpi_def_dim(ncid, "tot_num_g_cell_faces", tot_quants[NUM_CELL_FACES],
          &dimids[4]); ERR;
  err = ncmpi_def_dim(ncid, "tot_num_g_face_verts", tot_quants[NUM_FACE_VERTS],
          &dimids[5]); ERR;
  err = ncmpi_def_dim(ncid, "num_g_orig_particles", tot_quants[NUM_ORIG_PARTS],
          &dimids[6]); ERR;
  err = ncmpi_def_dim(ncid, "num_g_neighbors", tot_quants[NUM_NEIGHBORS],
          &dimids[7]); ERR;

  /* define variables */
  err = ncmpi_def_var(ncid, "num_verts", NC_INT, 1, &dimids[0],
          &varids[0]); ERR;
  err = ncmpi_def_var(ncid, "num_complete_cells", NC_INT, 1, &dimids[0],
          &varids[1]); ERR;
  err = ncmpi_def_var(ncid, "tot_num_cell_faces", NC_INT, 1, &dimids[0],
          &varids[2]); ERR;
  err = ncmpi_def_var(ncid, "tot_num_face_verts", NC_INT, 1, &dimids[0],
          &varids[3]); ERR;
  err = ncmpi_def_var(ncid, "num_orig_particles", NC_INT, 1, &dimids[0],
          &varids[4]); ERR;

  /* block offsets */
  err = ncmpi_def_var(ncid, "block_off_num_verts", NC_INT64, 1, &dimids[0],
          &varids[5]); ERR;
  err = ncmpi_def_var(ncid, "block_off_num_complete_cells", NC_INT64, 1,
          &dimids[0], &varids[6]); ERR;
  err = ncmpi_def_var(ncid, "block_off_tot_num_cell_faces", NC_INT64, 1,
          &dimids[0], &varids[7]); ERR;
  err = ncmpi_def_var(ncid, "block_off_tot_num_face_verts", NC_INT64, 1,
          &dimids[0], &varids[8]); ERR;
  err = ncmpi_def_var(ncid, "block_off_num_orig_particles", NC_INT64, 1,
          &dimids[0], &varids[9]); ERR;

  dimids_2D[0] = dimids[0];
  dimids_2D[1] = dimids[1];
  err = ncmpi_def_var(ncid, "mins", NC_FLOAT, 2, dimids_2D, &varids[11]); ERR;
  err = ncmpi_def_var(ncid, "maxs", NC_FLOAT, 2, dimids_2D, &varids[12]); ERR;

  dimids_2D[0] = dimids[2];
  dimids_2D[1] = dimids[1];
  err = ncmpi_def_var(ncid, "save_verts", NC_FLOAT, 2, dimids_2D,
          &varids[13]); ERR;
  dimids_2D[0] = dimids[6];
  dimids_2D[1] = dimids[1];
  err = ncmpi_def_var(ncid, "sites", NC_FLOAT, 2, dimids_2D,
          &varids[14]); ERR;
  err = ncmpi_def_var(ncid, "complete_cells", NC_INT, 1, &dimids[3],
          &varids[15]); ERR;
  err = ncmpi_def_var(ncid, "areas", NC_FLOAT, 1, &dimids[3],
          &varids[16]); ERR;
  err = ncmpi_def_var(ncid, "vols", NC_FLOAT, 1, &dimids[3], &varids[17]); ERR;
  err = ncmpi_def_var(ncid, "num_cell_faces", NC_INT, 1, &dimids[3],
          &varids[18]); ERR;
  err = ncmpi_def_var(ncid, "num_face_verts", NC_INT, 1, &dimids[4],
          &varids[19]); ERR;
  err = ncmpi_def_var(ncid, "face_verts", NC_INT, 1, &dimids[5],
          &varids[20]); ERR;
  err = ncmpi_def_var(ncid, "neighbors", NC_INT, 1, &dimids[7],
          &varids[21]); ERR;
  err = ncmpi_def_var(ncid, "g_block_ids", NC_INT, 1, &dimids[0],
          &varids[22]); ERR;

  /* exit define mode */
  err = ncmpi_enddef(ncid); ERR;

  /* write all variables.
     to improve: we can try nonblocking I/O to aggregate small requests */

  for (b = 0; b < nblocks; b++) {

    struct vblock_t *v = &vblocks[b];

    /* quantities */
    start[0] = block_ofsts[NUM_BLOCKS];
    count[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varids[0], start, count,
         &v->num_verts); ERR;
    err = ncmpi_put_vara_int_all(ncid, varids[1], start, count,
         &v->num_complete_cells); ERR;
    err = ncmpi_put_vara_int_all(ncid, varids[2], start, count,
         &v->tot_num_cell_faces); ERR;
    err = ncmpi_put_vara_int_all(ncid, varids[3], start, count,
         &v->tot_num_face_verts); ERR;
    err = ncmpi_put_vara_int_all(ncid, varids[4], start, count,
         &v->num_orig_particles); ERR;

    /* block offsets */
    err = ncmpi_put_vara_longlong_all(ncid, varids[5], start, count,
              &block_ofsts[NUM_VERTS]); ERR;
    err = ncmpi_put_vara_longlong_all(ncid, varids[6], start, count,
              &block_ofsts[NUM_COMP_CELLS]); ERR;
    err = ncmpi_put_vara_longlong_all(ncid, varids[7], start, count,
              &block_ofsts[NUM_CELL_FACES]); ERR;
    err = ncmpi_put_vara_longlong_all(ncid, varids[8], start, count,
              &block_ofsts[NUM_FACE_VERTS]); ERR;
    err = ncmpi_put_vara_longlong_all(ncid, varids[9], start, count,
              &block_ofsts[NUM_ORIG_PARTS]); ERR;

    /* block bounds */
    start[0] = block_ofsts[NUM_BLOCKS];
    count[0] = 1;
    start[1] = 0;
    count[1] = 3;
    err = ncmpi_put_vara_float_all(ncid, varids[11], start, count,
           v->mins); ERR;
    err = ncmpi_put_vara_float_all(ncid, varids[12], start, count,
           v->maxs); ERR;

    /* save_verts */
    start[0] = block_ofsts[NUM_VERTS];
    start[1] = 0;
    count[0] = v->num_verts;
    count[1] = 3;
    err = ncmpi_put_vara_float_all(ncid, varids[13], start, count,
           v->save_verts); ERR;

    /* sites */
    start[0] = block_ofsts[NUM_ORIG_PARTS];
    start[1] = 0;
    count[0] = v->num_orig_particles;
    count[1] = 3;
    err = ncmpi_put_vara_float_all(ncid, varids[14], start, count,
           v->sites); ERR;

    /* complete cells */
    start[0] = block_ofsts[NUM_COMP_CELLS];
    count[0] = v->num_complete_cells;
    err = ncmpi_put_vara_int_all(ncid, varids[15], start, count,
         v->complete_cells); ERR;

    /* areas */
    start[0] = block_ofsts[NUM_COMP_CELLS];
    count[0] = v->num_complete_cells;
    err = ncmpi_put_vara_float_all(ncid, varids[16], start, count,
           v->areas); ERR;

    /* volumes */
    start[0] = block_ofsts[NUM_COMP_CELLS];
    count[0] = v->num_complete_cells;
    err = ncmpi_put_vara_float_all(ncid, varids[17], start, count,
           v->vols); ERR;

    /* num_cell_faces */
    start[0] = block_ofsts[NUM_COMP_CELLS];
    count[0] = v->num_complete_cells;
    err = ncmpi_put_vara_int_all(ncid, varids[18], start, count,
         v->num_cell_faces); ERR;

    /* num_face_verts */
    start[0] = block_ofsts[NUM_CELL_FACES];
    count[0] = v->tot_num_cell_faces;
    err = ncmpi_put_vara_int_all(ncid, varids[19], start, count,
         v->num_face_verts); ERR;

    /* face verts */
    start[0] = block_ofsts[NUM_FACE_VERTS];
    count[0] = v->tot_num_face_verts;
    err = ncmpi_put_vara_int_all(ncid, varids[20], start, count,
         v->face_verts); ERR;

    /* neighbors */
    int *neighbors = (int*)malloc(DIY_Num_neighbors(0, b) * sizeof(int));
    int num_neighbors = DIY_Get_neighbors(0, b, neighbors);
    start[0] = block_ofsts[NUM_NEIGHBORS];
    count[0] = num_neighbors;
    err = ncmpi_put_vara_int_all(ncid, varids[21], start, count, neighbors);
    ERR;

    /* gids */
    int gid = DIY_Gid(0, b);
    start[0] = block_ofsts[NUM_BLOCKS];
    count[0] = 1;
    err = ncmpi_put_vara_int_all(ncid, varids[22], start, count,
         &gid); ERR;

    /* update block offsets */
    block_ofsts[NUM_VERTS] += v->num_verts;
    block_ofsts[NUM_COMP_CELLS] += v->num_complete_cells;
    block_ofsts[NUM_CELL_FACES] += v->tot_num_cell_faces;
    block_ofsts[NUM_FACE_VERTS] += v->tot_num_face_verts;
    block_ofsts[NUM_ORIG_PARTS] += v->num_orig_particles;
    block_ofsts[NUM_NEIGHBORS] += num_neighbors;
    block_ofsts[NUM_BLOCKS]++;

    /* debug */
/*     fprintf(stderr, "gid = %d num_verts = %d num_complete_cells = %d " */
/* 	    "tot_num_cell_faces = %d tot_num_face_verts = %d " */
/* 	    "num_orig_particles = %d\n", */
/* 	    gid, v->num_verts, v->num_complete_cells, v->tot_num_cell_faces, */
/* 	    v->tot_num_face_verts, v->num_orig_particles); */

  }

  err = ncmpi_close(ncid); ERR;
#endif

}
/*--------------------------------------------------------------------------*/
/*
  reads input in pnetcdf format

  nblocks: (output) local number of blocks
  tot_blocks: (output) total number of blocks
  vblocks: (output) pointer to array of vblocks
  in_file: input file name
  comm: MPI communicator
  gids: (output) gids of local blocks (allocated by this function)
  num_neighbors: (output) number of neighbors for each local block
   (allocated by this function)
  neighbors: (output) gids of neighbors of each local block
   (allocated by this function)

  side effects: allocates vblocks, gids, num_neighbors, neighbors

*/
void pnetcdf_read(int *nblocks, int *tot_blocks, struct vblock_t ***vblocks,
      char *in_file, MPI_Comm comm, int *gids, int *num_neighbors,
      int **neighbors) {

#ifdef USEPNETCDF
  int err;
  int ncid, varids[23], dimids[8];
  MPI_Offset start[2], count[2];
  nc_type type;
  int ndims, natts;
  int dims[2];
  int rank, groupsize; /* MPI usual */

  /* open file for reading */
  err = ncmpi_open(comm, in_file, NC_NOWRITE, MPI_INFO_NULL, &ncid); ERR;

  err = ncmpi_inq_varid(ncid, "block_off_num_verts", &varids[5]); ERR;
  err = ncmpi_inq_varid(ncid, "block_off_num_complete_cells", &varids[6]); ERR;
  err = ncmpi_inq_varid(ncid, "block_off_tot_num_cell_faces", &varids[7]); ERR;
  err = ncmpi_inq_varid(ncid, "block_off_tot_num_face_verts", &varids[8]); ERR;
  err = ncmpi_inq_varid(ncid, "block_off_num_orig_particles", &varids[9]); ERR;

  /* get number of blocks */
  MPI_Offset num_g_blocks; /* 64 bit version of tot_blcoks */
  err = ncmpi_inq_dimlen(ncid, dimids[0], &num_g_blocks); ERR;
  *tot_blocks = num_g_blocks;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);
  int start_block_ofst =  rank * (*tot_blocks / groupsize);
  *nblocks = (rank < groupsize - 1 ? (*tot_blocks / groupsize) :
        *tot_blocks - (rank * *tot_blocks / groupsize));

  /* block offsets */
  int64_t *block_ofsts = (int64_t*)malloc(*tot_blocks * sizeof(int64_t));
  *vblocks = (struct vblock_t**)malloc(*nblocks * sizeof(struct vblock_t*));

  /* read all blocks */
  gids = (int *)malloc(*nblocks * sizeof(int));
  num_neighbors = (int *)malloc(*nblocks * sizeof(int));
  neighbors = (int **)malloc(*nblocks * sizeof(int *));
  int b;
  for (b = 0; b < *nblocks; b++) {

    struct vblock_t* v = (struct vblock_t*)malloc(sizeof(struct vblock_t));

    /* quantities */
    start[0] = start_block_ofst + b;
    count[0] = 1;
    err = ncmpi_inq_varid(ncid, "num_verts", &varids[0]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[0], start, count,
         &(v->num_verts)); ERR;
    err = ncmpi_inq_varid(ncid, "num_complete_cells", &varids[1]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[1], start, count,
         &(v->num_complete_cells)); ERR;
    err = ncmpi_inq_varid(ncid, "tot_num_cell_faces", &varids[2]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[2], start, count,
         &(v->tot_num_cell_faces)); ERR;
    err = ncmpi_inq_varid(ncid, "tot_num_face_verts", &varids[3]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[3], start, count,
         &(v->tot_num_face_verts)); ERR;
    err = ncmpi_inq_varid(ncid, "num_orig_particles", &varids[4]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[4], start, count,
         &(v->num_orig_particles)); ERR;
    err = ncmpi_inq_varid(ncid, "neighbors", &varids[21]); ERR;
    err = ncmpi_inq_var(ncid, varids[21], 0, &type, &ndims,
      dimids, &natts);

    /* block bounds */
    start[0] = start_block_ofst + b;
    start[1] = 0;
    count[0] = 1;
    count[1] = 3;
    err = ncmpi_inq_varid(ncid, "mins", &varids[11]); ERR;
    err = ncmpi_get_vara_float_all(ncid, varids[11], start, count,
           v->mins); ERR;
    err = ncmpi_inq_varid(ncid, "maxs", &varids[12]); ERR;
    err = ncmpi_get_vara_float_all(ncid, varids[12], start, count,
           v->maxs); ERR;

    /* save_verts */
    start[0] = 0;
    count[0] = *tot_blocks;
    err = ncmpi_get_vara_longlong_all(ncid, varids[5], start, count,
              block_ofsts); ERR;
    v->save_verts = (float *)malloc(v->num_verts * 3 * sizeof(float));
    start[0] = block_ofsts[start_block_ofst + b];
    start[1] = 0;
    count[0] = v->num_verts;
    count[1] = 3;
    err = ncmpi_inq_varid(ncid, "save_verts", &varids[13]); ERR;
    err = ncmpi_get_vara_float_all(ncid, varids[13], start, count,
           v->save_verts); ERR;

    /* sites */
    start[0] = 0;
    count[0] = *tot_blocks;
    err = ncmpi_get_vara_longlong_all(ncid, varids[9], start, count,
              block_ofsts); ERR;
    v->sites = (float *)malloc(v->num_orig_particles * 3 * sizeof(float));
    start[0] = block_ofsts[start_block_ofst + b];
    start[1] = 0;
    count[0] = v->num_orig_particles;
    count[1] = 3;
    err = ncmpi_inq_varid(ncid, "sites", &varids[14]); ERR;
    err = ncmpi_get_vara_float_all(ncid, varids[14], start, count,
           v->sites); ERR;

    /* complete cells */
    start[0] = 0;
    count[0] = *tot_blocks;
    err = ncmpi_get_vara_longlong_all(ncid, varids[6], start, count,
              block_ofsts); ERR;
    v->complete_cells = (int *)malloc(v->num_complete_cells * sizeof(int));
    start[0] = block_ofsts[start_block_ofst + b];
    count[0] = v->num_complete_cells;
    err = ncmpi_inq_varid(ncid, "complete_cells", &varids[15]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[15], start, count,
         v->complete_cells); ERR;

    /* areas, uses same block offsets as complete cells */
    v->areas = (float *)malloc(v->num_complete_cells * sizeof(float));
    start[0] = block_ofsts[start_block_ofst + b];
    count[0] = v->num_complete_cells;
    err = ncmpi_inq_varid(ncid, "areas", &varids[16]); ERR;
    err = ncmpi_get_vara_float_all(ncid, varids[16], start, count,
           v->areas); ERR;

    /* volumes, uses same block offsets as complete cells */
    v->vols = (float *)malloc(v->num_complete_cells * sizeof(float));
    start[0] = block_ofsts[start_block_ofst + b];
    count[0] = v->num_complete_cells;
    err = ncmpi_inq_varid(ncid, "vols", &varids[17]); ERR;
    err = ncmpi_get_vara_float_all(ncid, varids[17], start, count,
           v->vols); ERR;

    /* num_cell_faces, uses same block offsets as complete cells */
    v->num_cell_faces = (int *)malloc(v->num_complete_cells * sizeof(int));
    start[0] = block_ofsts[start_block_ofst + b];
    count[0] = v->num_complete_cells;
    err = ncmpi_inq_varid(ncid, "num_cell_faces", &varids[18]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[18], start, count,
         v->num_cell_faces); ERR;

    /* num_face_verts */
    start[0] = 0;
    count[0] = *tot_blocks;
    err = ncmpi_get_vara_longlong_all(ncid, varids[7], start, count,
              block_ofsts); ERR;
    v->num_face_verts = (int *)malloc(v->tot_num_cell_faces * sizeof(int));
    start[0] = block_ofsts[start_block_ofst + b];
    count[0] = v->tot_num_cell_faces;
    err = ncmpi_inq_varid(ncid, "num_face_verts", &varids[19]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[19], start, count,
         v->num_face_verts); ERR;

    /* face_verts */
    start[0] = 0;
    count[0] = *tot_blocks;
    err = ncmpi_get_vara_longlong_all(ncid, varids[8], start, count,
              block_ofsts); ERR;
    v->face_verts = (int *)malloc(v->tot_num_face_verts * sizeof(int));
    start[0] = block_ofsts[start_block_ofst + b];
    count[0] = v->tot_num_face_verts;
    err = ncmpi_inq_varid(ncid, "face_verts", &varids[20]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[20], start, count,
         v->face_verts); ERR;

    /* neighbors */
    MPI_Offset n; /* temporary 64-bit version of number of neighbors */
    err = ncmpi_inq_varid(ncid, "neighbors", &varids[21]); ERR;
    err = ncmpi_inq_var(ncid, varids[2], 0, &type, &ndims,
      dims, &natts); ERR;
    err = ncmpi_inq_dimlen(ncid, dims[0], &n); ERR;
    num_neighbors[b] = n;
    neighbors[b] = (int *)malloc(num_neighbors[b] * sizeof(int));
    start[0] = start_block_ofst + b;
    count[0] = num_neighbors[b];
    err = ncmpi_get_vara_int_all(ncid, varids[21], start, count,
         neighbors[b]); ERR;

    /* gids */
    start[0] = start_block_ofst + b;
    count[0] = 1;
    err = ncmpi_inq_varid(ncid, "g_block_ids", &varids[22]); ERR;
    err = ncmpi_get_vara_int_all(ncid, varids[22], start, count,
         &gids[b]); ERR;

    (*vblocks)[b] = v;

  }

  /* cleanup */
  err = ncmpi_close(ncid); ERR;
  free(block_ofsts);
#endif
}
/*--------------------------------------------------------------------------*/
/*
  creates DIY datatype for the subset of the voronoi block to write to disk

  vblock: voronoi block
  did: domain id (unused in this case)
  lid: local block number (unused in this case)
  dtype: pointer to datatype

  side effects: commits datatype

*/
void create_datatype(void* vblock, int did, int lid, DIY_Datatype *dtype) {

  did = did; /* quiet compiler warning */
  lid = lid;

  struct vblock_t *v = (struct vblock_t *)vblock;
  struct map_block_t map[] = {

    { DIY_FLOAT,  OFST, 3,
      offsetof(struct vblock_t, mins)                 },
    { DIY_FLOAT, ADDR, v->num_verts * 3,
      DIY_Addr(v->save_verts)                              },
    { DIY_FLOAT, ADDR, v->num_orig_particles * 3,
      DIY_Addr(v->sites)                              },
    { DIY_INT,    ADDR, v->num_complete_cells,
      DIY_Addr(v->complete_cells)                     },
    { DIY_FLOAT, ADDR, v->num_complete_cells,
      DIY_Addr(v->areas)                              },
    { DIY_FLOAT, ADDR, v->num_complete_cells,
      DIY_Addr(v->vols)                               },
    { DIY_INT,    ADDR, v->num_complete_cells,
      DIY_Addr(v->num_cell_faces)                     },
    { DIY_INT,    ADDR, v->tot_num_cell_faces,
      DIY_Addr(v->num_face_verts)                     },
    { DIY_INT,    ADDR, v->tot_num_face_verts,
      DIY_Addr(v->face_verts)                         },
    { DIY_FLOAT,  OFST, 3,
      offsetof(struct vblock_t, maxs)                 },

  };

  DIY_Create_struct_datatype(DIY_Addr(vblock), 10, map, dtype);

#ifdef DATATYPE_PROFILE
  /* print datatype data for John Jenkins */
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  fprintf(stderr, "Rank %d creating MPI struct datatype:\n"
    "base type = MPI_INT,   quantity = 7, address = %p\n"
    "base type = MPI_FLOAT, quantity = 3, address = %p\n"
    "base type = MPI_FLOAT, quantity = %d, address = %p\n"
    "base type = MPI_FLOAT, quantity = %d, address = %p\n"
    "base type = MPI_INT,   quantity = %d, address = %p\n"
    "base type = MPI_FLOAT, quantity = %d, address = %p\n"
    "base type = MPI_FLOAT, quantity = %d, address = %p\n"
    "base type = MPI_INT,   quantity = %d, address = %p\n"
    "base type = MPI_INT,   quantity = %d, address = %p\n"
    "base type = MPI_INT,   quantity = %d, address = %p\n"
    "base type = MPI_FLOAT, quantity = 3, address =  %p\n\n\n",
    rank,
    hdrs[lid],
    vblock + offsetof(struct vblock_t, mins),
    v->num_verts * 3, v->save_verts,
    v->num_orig_particles * 3, v->sites,
    v->num_complete_cells, v->complete_cells,
    v->num_complete_cells, v->areas,
    v->num_complete_cells, v->vols,
    v->num_complete_cells, v->num_cell_faces,
    v->tot_num_cell_faces, v->num_face_verts,
    v->tot_num_face_verts, v->face_verts,
    vblock + offsetof(struct vblock_t, maxs));
#endif

}
/*--------------------------------------------------------------------------*/

