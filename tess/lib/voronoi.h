/*---------------------------------------------------------------------------
 *
 * voronoi data structures
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
#ifndef _VORONOI_H
#define _VORONOI_H

/* header */
/* #define HDR_ELEMS 7 /\* number of header elements used *\/ */
#define NUM_VERTS 0 /* indices of header elements */
#define NUM_CELLS 1
#define TOT_NUM_CELL_VERTS 2
#define NUM_COMPLETE_CELLS 3
#define TOT_NUM_CELL_FACES 4
#define TOT_NUM_FACE_VERTS 5
#define NUM_ORIG_PARTICLES 6

#define AVG_CELL_FACES 32 /* average number of faces per cell */
#define AVG_FACE_VERTS 16 /* average number of vertices per face */
#define MAX_HIST_BINS 256 /* maximum number of bins in cell volume histogram */

/* timing info */
#define MAX_TIMES 8
#define EXCH_TIME 0
#define CELL_TIME 1
#define HULL_TIME 2
#define OUT_TIME  3

/* voronoi diagram for one block */
struct vblock_t {

  float mins[3]; /* block minimum corner */

  /* raw voronoi tesselation */
  int num_verts; /* number of vertices, including "infinite" vertex */
  int num_cells; /* number of voronoi cells, same as number of input sites */
  double *verts; /* double version of vertices x0, y0, z0, x1, y, z1, ...
		    infinite vertex is x0, y0, z0, for in-core computations */
  float *save_verts; /* float version of verts for permananent storage */
  int *num_cell_verts; /* number of vertices in each cell */
  int tot_num_cell_verts; /* sum of all num_cell_verts */
  int *cells; /* vertices in cells, in order of cells */
  float *sites; /* original input points = sites of voronoi cells
		    in same order as cells */

  /* temporary complete cells */
  int temp_num_complete_cells; /* number of complete cells */
  int *temp_complete_cells; /* sorted list of complete cells from 
			       raw tesselation */

  /* final complete cells */
  int num_complete_cells; /* number of complete cells */
  int *complete_cells; /* sorted list of complete cells from raw tesselation */
  float *areas; /* surface areas of complete cells */
  float *vols; /* volumes of complete cells */
  int tot_num_cell_faces; /* total number of faces in complete cells */
  int *num_cell_faces; /* number of faces in complete cells, in order of
			  complete cells */
  int *num_face_verts; /* number of vertices in each face of complete cells, 
			  in order of complete cells and faces in them 
			  num_face_verts[comp_cell][face] */
  int tot_num_face_verts; /* total number of vertices in faces of 
			     complete cells */  
  int *face_verts; /* vertices in each face of complete cells, 
		      in order of complete cells and faces of in them 
		      face_verts[comp_cell][face][vertex] */

  /* additional info */
  int num_orig_particles; /* number of original particles in this block
			     before any neighbor exhcange */
  float maxs[3]; /* block maximum corner */

};

/* temporary convex hull output for one cell */
struct cblock_t {

  float area; /* surface area */
  float vol; /* volume */
  int num_cell_faces; /* number of faces */
  int *num_face_verts; /* number of vertices in each face */
  int **face_verts; /* vertices in each face */

};

/* statistical summary */
struct stats_t {
  int tot_cells; /* total number of complete voronoi cells found */
  int tot_faces; /* total number of faces in all complete voronoi cells */
  int tot_verts; /* total number of vertices in all complete voronoi cells */
  float avg_cell_verts; /* average number of vertices per cell */
  float avg_cell_faces; /* average number of faces per cell */
  float avg_face_verts; /* average number of vertices per face */
  float min_cell_vol; /* minimum cell volume */
  float max_cell_vol; /* maximum cell volume */
  float avg_cell_vol; /* average cell volume */
  float min_cell_dense; /* minimum cell density */
  float max_cell_dense; /* maximum cell density */
  float avg_cell_dense; /* average cell density */
  float min_cell_time; /* minimum voronoi cell time */
  float max_cell_time; /* maximum voronoi cell time */
  float min_hull_time; /* minimum convex hull time */
  float max_hull_time; /* maximum convex hull time */
  int num_vol_bins; /* number of bins in cell volume histogram */
  int num_dense_bins; /* number of bins in cell density histogram */
  int vol_hist[MAX_HIST_BINS]; /* cell volume histogram */
  int dense_hist[MAX_HIST_BINS]; /* cell density histogram */
};

#endif
