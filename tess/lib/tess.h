/*---------------------------------------------------------------------------
 *
 * parallel tesselation
 *
 * Tom Peterka
 * Argonne National Laboratory
 * 9700 S. Cass Ave.
 * Argonne, IL 60439
 * tpeterka@mcs.anl.gov
 *
 * (C) 2011 by Argonne National Laboratory.
 * See COPYRIGHT in top-level directory.
 *
--------------------------------------------------------------------------*/
#ifndef _TESS_H
#define _TESS_H

#define qh_QHimport
#include "qhull_a.h"

#include "voronoi.h"

/* public */

#ifdef __cplusplus
extern "C"
#endif
void tess_test(int tot_blocks, int *data_size, float jitter, float cell_size,
	       float ghost_factor, float minvol, float maxvol, 
	       double *all_times);

#ifdef __cplusplus
extern "C"
#endif
void tess_init(int num_blocks, int *gids, 
	       struct bb_t *bounds, struct gb_t **neighbors, 
	       int *num_neighbors, float cell_dia, float ghost_mult, 
	       float *global_mins, float *global_maxs, 
	       int wrap, float minvol, float maxvol, MPI_Comm mpi_comm,
	       double *times);
#ifdef __cplusplus
extern "C"
#endif
void tess_finalize();

#ifdef __cplusplus
extern "C"
#endif
void tess(float **particles, int *num_particles, char *out_file);

/* private */

void voronoi(int nblocks, float **particles, int *num_particles, 
	     float cell_size, float ghost_factor,
	     double *times, char *out_file);
int gen_particles(int lid, float **particles, float jitter);
void gen_voronoi_output(facetT *facetlist, struct vblock_t *vblock);
void *create_datatype(void *vblock, int lid, MPI_Datatype *dtype);
void test_datatype(struct vblock_t *vblock, int *hdr, 
		   void* (*type_func)(void*, MPI_Datatype*));
int gen_convex_output(facetT *facetlist, struct cblock_t *cblock);
void convex_to_voronoi(struct cblock_t *cblock, struct vblock_t *vblock,
		       int *vmap, int cell);
void complete_cells(struct vblock_t *vblock, int lid, float ghost);
void prep_vertices(struct vblock_t *vblock, int comp_cell, double *vertices,
		   int *vmap);
void create_blocks(int num_blocks, struct vblock_t **vblocks, int ***hdrs);
void orig_cells(int nblocks, struct vblock_t *vblocks, int dim,
		int *num_particles, int *num_orig_particles, 
		float **particles, float ghost);
void cell_hulls(int nblocks, struct vblock_t *vblocks, int dim);
void neighbor_particles(int nblocks, float **particles, int *num_particles, 
			float ghost);
void write_out(int nblocks, struct vblock_t *vblocks, int **hdrs,
	       char *out_file);
MPI_Datatype* recv_type(int *cts);
MPI_Datatype* send_type(int *cts, char** pts);
void del_recv_type(MPI_Datatype *dtype);
void del_send_type(MPI_Datatype *dtype);
void destroy_blocks(int num_blocks, struct vblock_t *vblocks, int **hdrs);
void collect_stats(int nblocks, struct vblock_t *vblocks, double *times);
void aggregate_stats(int nblocks, struct vblock_t *vblocks, 
		     struct stats_t *loc_stats);
void average(void *in, void *inout, int *len, MPI_Datatype *type);
void histogram(void *in, void *inout, int *len, MPI_Datatype *type);
void print_block(struct vblock_t *vblock, int gid);
void print_particles(float *particles, int num_particles, int gid);
void prep_out(int nblocks, struct vblock_t *vblocks, int **hdrs);
void transform_particle(char *p, unsigned char wrap_dir);

#endif
