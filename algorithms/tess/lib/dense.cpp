//---------------------------------------------------------------------------
//
// density field regular grid computation from voronoi tessellation
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// (C) 2013 by Argonne National Laboratory.
// See COPYRIGHT in top-level directory.
//
//--------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <vector>
#include "voronoi.h"
#include <math.h>
#include "mpi.h"
#include "diy.h"

using namespace std;

// voronoi blocks
vblock_t **vblocks;

// grid point
struct grid_pt_t {
  int idx[3]; // global grid point index
  float dense; // density
};

// debug
float min_cell_vol = 1000.0; // min and max cell volume
float max_cell_vol = 0.0;
float min_sampled_cell_vol = 1000.0; // min and max cell volume w/ sample pts
float max_sampled_cell_vol = 0.0;
float min_dense = 1000.0; // min and max density
float max_dense = 0.0;
int max_idx = 0; // maximum idx written

// globals
static float data_mins[3], data_maxs[3]; // data global extents
static float grid_step_size[3]; // physical size of one grid space
static int glo_num_idx[3]; // global grid number of points
static int swap_bytes; // swap bytes flag
static char infile[256]; // input file name
static float mass; // particle mass
static char outfile[] = "dense.raw"; // output file name
static float eps = 0.0001; // epsilon for floating point values to be equal

// function prototypes
void ParseArgs(int argc, char **argv);
void BlockGridParams(int lid, int *block_min_idx, int *block_max_idx,
		     int *block_num_idx);
void IterateCells(int block, int *block_num_idx, float **density);

// DEPRECATED
// but not commented out until sure that weighted scheme not needed
void IterateCellsWeighted(int block, int *block_min_idx, int *block_max_idx,
			  int *block_num_idx, float **density);

void *CreateReadType(int did, int lid, int *hdr, DIY_Datatype *dtype);
void CellBounds(vblock_t *vblock, int &cell, int &face, int &vert,
		float *cell_min, float *cell_max, float *centroid);
void CellGridBounds(float *cell_mins, float *cell_maxs, 
		    int *block_min_idx, int *block_max_idx,
		    int *num_cell_grid_pts, 
		    int *cell_min_grid_idx, int *cell_max_grid_idx);
void CellGridParams(int *cell_grid_pts, int *cell_min_grid_idx,
		    float *cell_min_grid_pos, float *cell_min, 
		    float *cell_max);
int CellGridPointParams(int *grid_idx, int *block_grid_idx, float *grid_pos,
			int xi, int yi, int zi, int *cell_min_grid_idx,
			int *block_min_grid_idx, float *cell_min_grid_pos);
void DistributeScalar(float *pt, float scalar, float *win_min, float *win_max, 
		      int *block_min_idx, int *block_max_idx,
		      vector <int> &grid_idxs, vector <float> &grid_scalars);
void ComputeNormal(float *verts, int num_verts, float *normal);
float PtPlaneDist(float *pt, float *verts, int num_verts);
bool PtInCell(float *pt, vblock_t *vblock, int &cell, int &face, int &vert);
void Global2LocalIdx(int *global_idx, int *local_idx, int block);
int CellGridPts(int block, vblock_t *vblock, int &cell, int &face, int &vert,
		float *cell_mins, float *cell_maxs, 
		vector <grid_pt_t> &grid_pts);

// DEPRECATED
// void DistributeDensity(float *pt, float dense, float win_frac, 
// 		       vector <int> &grid_idxs, vector <float> &grid_denses);

void ItemDtype(DIY_Datatype *dtype);
void WriteGrid(float **density, MPI_Comm comm, int nblocks, 
	       int mblocks);
void handle_error(int errcode, char *str, MPI_Comm comm);
int index(int *block_grid_idx, int *block_num_idx);
bool in_block(int *idx, int *min_idx, int *max_idx);
void phys_pos(int *grid_idx, float *grid_pos);

//--------------------------------------------------------------------------

int main(int argc, char** argv) {

  int dim = 3;
  int tot_blocks; // global number of blocks
  int nblocks; // my local number of blocks
  int num_threads = 1; // number of threads DIY can use
  int rank, groupsize; // MPI usual
  int did; // domain id
  MPI_Comm comm = MPI_COMM_WORLD; // MPI communicator

  // local block grid parameters
  int block_min_idx[3]; // global grid index of block minimum grid point
  int block_max_idx[3]; // global grid index of block maximum grid point
  int block_num_idx[3]; // number of grid points in local block

  // parse args
  ParseArgs(argc, argv);

  // init
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &groupsize);
  // data_pts for initializing DIY doesn't really matter in this case
  // could be anything, just make it based on unit spacing by defualt
  int data_pts[3] = {data_maxs[0] - data_mins[0] + 1,
		     data_maxs[1] - data_mins[1] + 1,
		     data_maxs[2] - data_mins[2] + 1}; 
  DIY_Init(dim, data_pts, num_threads, comm);
  did = DIY_Read_decomposed(infile, swap_bytes, &tot_blocks, &nblocks);
  int maxblocks = (int)(ceil(tot_blocks / (float)groupsize));

  // read voronoi blocks
  int **hdrs = new int*[nblocks];
  for (int i = 0; i < nblocks; i++)
    hdrs[i] = new int[DIY_MAX_HDR_ELEMENTS];
  DIY_Read_open_all(did, infile, swap_bytes, 0, &tot_blocks, &nblocks);
  DIY_Read_blocks_all(did, (void ***)(&vblocks), hdrs, &CreateReadType);
  DIY_Read_close_all(did);

  // debug
  for (int i = 0; i < nblocks; i++) {
    bb_t bb; // block bounds
    DIY_Block_bounds(did, i, &bb);
    fprintf(stderr, "gid[%d] = %d DIY bounds: min = [%.3f %.3f %.3f] "
	    "max[%.3f %.3f %.3f] "
	    "Voronoi bounds: min = [%.3f %.3f %.3f] "
	    "max[%.3f %.3f %.3f]\n", 
	    i, DIY_Gid(did, i), bb.min[0], bb.min[1], bb.min[2],
	    bb.max[0], bb.max[1], bb.max[2],
	    vblocks[i]->mins[0], vblocks[i]->mins[1], vblocks[i]->mins[2],
	    vblocks[i]->maxs[0], vblocks[i]->maxs[1], vblocks[i]->maxs[2]);
  }

  // allocate density field
  float *density[nblocks];
  for (int block = 0; block < nblocks; block++) {
    BlockGridParams(block, block_min_idx, block_max_idx, block_num_idx);
    int npts = block_num_idx[0] * block_num_idx[1] * block_num_idx[2];
    density[block] = new float[npts];
    memset(density[block], 0 , npts * sizeof(float));

    // debug
    fprintf(stderr, "The grid in block %d extends from [%d %d %d] to "
	    "[%d %d %d] and is [%d %d %d] in size.\n",
	    block, block_min_idx[0], block_min_idx[1], block_min_idx[2],
	    block_max_idx[0], block_max_idx[1], block_max_idx[2],
	    block_num_idx[0], block_num_idx[1], block_num_idx[2]);

  }

  // sample the density

  for (int block = 0; block < nblocks; block++) { // blocks

    // get local block grid parameters
    BlockGridParams(block, block_min_idx, block_max_idx, block_num_idx);

    // iterate over cells, distributing density onto grid points
    IterateCells(block, block_num_idx, density);

  }

  // the following is turned off until the domain decomposition properly
  // matches the file (todo)

#if 0

  // received items from neighbors
  void ***items = new void**[nblocks]; // received items
  int *num_items = new int[nblocks]; // number of received items in each block

  // exchange neighbors
  DIY_Exchange_neighbors(did, items, num_items, 1.0, &ItemDtype);

  // save received items
  for (int block = 0; block < nblocks; block++) {

    BlockGridParams(block, block_min_idx, block_max_idx, block_num_idx);

    // debug
//     fprintf(stderr, "%d items received by global block %d:\n", 
// 	    num_items[block], DIY_Gid(did, block));
//     fprintf(stderr, "block %d block_min_idx = [%d %d %d]\n", block,
// 	    block_min_idx[0], block_min_idx[1], block_min_idx[2]);

    for (int j = 0; j < num_items[block]; j++) {

      grid_pt_t *grid_pt = (grid_pt_t*)(items[block][j]);

      // debug
//       fprintf(stderr, "density[%d][%d][%d] = %.3f\n",
// 	      grid_pt->idx[0], grid_pt->idx[1], grid_pt->idx[2], 
// 	      grid_pt->dense);

      int block_grid_idx[3]; // indices in local block array
      block_grid_idx[0] = grid_pt->idx[0] - block_min_idx[0];
      block_grid_idx[1] = grid_pt->idx[1] - block_min_idx[1];
      block_grid_idx[2] = grid_pt->idx[2] - block_min_idx[2];

      // debug
//       fprintf(stderr, "saving [%d][%d][%d]\n", xi, yi, zi);

      int idx = index(block_grid_idx, block_num_idx);
      density[block][idx] = grid_pt->dense;

      // debug
      if (grid_pt->dense < min_dense)
	min_dense = grid_pt->dense;
      if (grid_pt->dense > max_dense)
	max_dense = grid_pt->dense;
      if (idx > max_idx)
	max_idx = idx;

    }

  }

#endif

  // debug
  fprintf(stderr, "min_cell_vol = %.3e max_cell_voll = %.3e\n"
	  "min_sampled_cell_vol = %.3e max_sampled_cell_voll = %.3e\n"
	  "min_dense = %.3e max_dense = %.3e\n",
	  min_cell_vol, max_cell_vol, min_sampled_cell_vol, 
	  max_sampled_cell_vol, min_dense, max_dense);

  // write file
  WriteGrid(density, comm, nblocks, maxblocks);

  // cleanup
  for (int i = 0; i < nblocks; i++)
    delete[] density[i];

  // the following is turned off until the domain decomposition properly
  // matches the file (todo)

#if 0

  DIY_Flush_neighbors(did, items, num_items, &ItemDtype);
  delete[] num_items;
  delete[] items;

#endif

  DIY_Finalize();
  MPI_Finalize();

}
//--------------------------------------------------------------------------
//
// iterate over cells and assign single density to grid point
//
// block: local block number
// block_num_idx: number of grid points in block (output) (x,y,z)
// density: density field
//
// side effects: writes density or sends to neighbors
//
void IterateCells(int block, int *block_num_idx, float **density) {

  int face = 0; // faces counter
  int vert = 0; // face vertices counter
  float cell_min[3], cell_max[3], cell_centroid[3]; // cell bounds
  float grid_pos[3]; // physical position of grid point

  // cells
  for (int cell = 0; cell < vblocks[block]->num_complete_cells; cell++) {

    // cell bounds
    CellBounds(vblocks[block], cell, face, vert, cell_min, cell_max,
	       cell_centroid);

    // debug
    if (vblocks[block]->vols[cell] < min_cell_vol)
      min_cell_vol = vblocks[block]->vols[cell];
    if (vblocks[block]->vols[cell] > max_cell_vol)
      max_cell_vol = vblocks[block]->vols[cell];

    // grid points covered by this cell
    vector <grid_pt_t> grid_pts;
    CellGridPts(block, vblocks[block], cell, face, vert, cell_min, cell_max, 
		grid_pts);

    // iterate over grid points covered by cell
    for (int i = 0; i < (int)grid_pts.size(); i++) {

      phys_pos(grid_pts[i].idx, grid_pos);

      // assign density to grid points in the block
      if (grid_pos[0] >= vblocks[block]->mins[0] &&
	  (grid_pos[0] < vblocks[block]->maxs[0]  ||
	   fabs(grid_pos[0] - data_maxs[0]) < eps) &&

	  grid_pos[1] >= vblocks[block]->mins[1] &&
	  (grid_pos[1] < vblocks[block]->maxs[1]  ||
	   fabs(grid_pos[1] - data_maxs[1]) < eps) &&

	  grid_pos[2] >= vblocks[block]->mins[2] &&
	  (grid_pos[2] < vblocks[block]->maxs[2]  ||
	   fabs(grid_pos[2] - data_maxs[2]) < eps) ) {

	// assign the density to the local block array
	int block_grid_idx[3]; // local block idx of grid point
	Global2LocalIdx(grid_pts[i].idx, block_grid_idx, block);
	int idx = index(block_grid_idx, block_num_idx);
	density[block][idx] = grid_pts[i].dense;

	// debug
	if (density[block][idx] < min_dense)
	  min_dense = density[block][idx];
	if (density[block][idx] > max_dense)
	  max_dense = density[block][idx];
	if (idx > max_idx)
	  max_idx = idx;

      }

      // the following is turned off until the domain decomposition properly
      // matches the file (todo)

#if 0

      // or send grid points to neighboring blocks
      else {

	// enqueue the point idx and density value and exchange
	DIY_Enqueue_item_points(did, block, (void *)&grid_pt, NULL,
				sizeof(grid_pt_t), grid_pos, 1, NULL);

	// debug
	if (grid_pt.dense < min_dense)
	  min_dense = grid_pt.dense;
	if (grid_pt.dense > max_dense)
	  max_dense = grid_pt.dense;

      }

#endif

    } // grid points covered by cell

  } // cells

}
//--------------------------------------------------------------------------
//
// DEPRECATED
// but not commented out until sure that weighted scheme not needed
//
// iterate over cells and assign average density to grid point
//
// block: local block number
// block_min_idx: global grid idx of block minimum grid point (x,y,z)
// block_max_idx: global grid idx of block maximum grid point (x,y,z)
// block_num_idx: number of grid points in block (x,y,z)
// density: density field
//
// side effects: writes density or sends to neighbors
//
void IterateCellsWeighted(int block, int *block_min_idx, int *block_max_idx,
			  int *block_num_idx, float **density) {

  int face = 0; // faces counter
  int vert = 0; // face vertices counter
  float cell_min[3], cell_max[3], cell_centroid[3]; // cell bounds

  // cells
  for (int cell = 0; cell < vblocks[block]->num_complete_cells; cell++) {

    // cell bounds
    CellBounds(vblocks[block], cell, face, vert, cell_min, cell_max,
	       cell_centroid);

    // debug
    if (vblocks[block]->vols[cell] < min_cell_vol) {
      min_cell_vol = vblocks[block]->vols[cell];
      min_sampled_cell_vol = vblocks[block]->vols[cell];
    }
    if (vblocks[block]->vols[cell] > max_cell_vol) {
      max_cell_vol = vblocks[block]->vols[cell];
      max_sampled_cell_vol = vblocks[block]->vols[cell];
    }

    // distribute density at cell centroid to neighboring grid points

    vector<int> grid_idxs; // grid idxs that get a fraction of the density
    vector<float> grid_denses; // density given to each grid_idx
    float dense = mass / vblocks[block]->vols[cell]; // total density of cell
    DistributeScalar(cell_centroid, dense, cell_min, cell_max, 
		     block_min_idx, block_max_idx, grid_idxs, 
		     grid_denses);
    for (int i = 0; i < (int)grid_idxs.size() / 3; i++) {

      // grid point is in this block
      if (in_block(&grid_idxs[3 * i], block_min_idx, block_max_idx)) {
	int block_idx[3]; // local block idx
	block_idx[0] = grid_idxs[3 * i] - block_min_idx[0];
	block_idx[1] = grid_idxs[3 * i + 1] - block_min_idx[1];
	block_idx[2] = grid_idxs[3 * i + 2] - block_min_idx[2];
	int idx = index(block_idx, block_num_idx);
	density[block][idx] += grid_denses[i];

      }

      // grid point is in another block
      // the following is turned off until the domain decomposition properly
      // matches the file (todo)

#if 0

      else {

	// enqueue the point idx and density value and exchange
	grid_pt_t grid_pt;
	grid_pt.idx[0] = grid_idxs[3 * i];
	grid_pt.idx[1] = grid_idxs[3 * i + 1];
	grid_pt.idx[2] = grid_idxs[3 * i + 2];
	grid_pt.dense = grid_denses[i];
	DIY_Enqueue_item_points(did, block, (void *)&grid_pt, NULL,
				sizeof(grid_pt_t), grid_pos, 1, NULL);

	// todo: don't forget to accumulate the received densities, 
	// individual samples, not accumlated totals are being sent

      }

#endif

    }

  } // cells

  // debug
  for (int i = 0; i < block_num_idx[0] * block_num_idx[1] * block_num_idx[2]; 
       i++) {
    if (density[block][i] < min_dense)
      min_dense = density[block][i];
    if (density[block][i] > max_dense)
      max_dense = density[block][i];
    if (i > max_idx)
      max_idx = i;
    if (i < 31)
      fprintf(stderr, "density[%d][%d] = %.2e\n", block, i, density[block][i]);
  }

}
//--------------------------------------------------------------------------
//
// allocates a new block and creates DIY datatype for it
//
// did: domain id
// lid: local block number
// hdr: block header
// dtype: pointer to datatype
//
// side effects: commits datatype, DIY will cleanup
//
// returns: address of new block
//
void *CreateReadType(int did, int lid, int *hdr, DIY_Datatype *dtype) {

  did = did; // quiet compiler warning
  lid = lid; 

  // allocate space

  struct vblock_t *v = new vblock_t;

  v->num_verts = hdr[NUM_VERTS];
  v->num_cells = hdr[NUM_CELLS];
  v->tot_num_cell_verts = hdr[TOT_NUM_CELL_VERTS];
  v->num_complete_cells = hdr[NUM_COMPLETE_CELLS];
  v->tot_num_cell_faces = hdr[TOT_NUM_CELL_FACES];
  v->tot_num_face_verts = hdr[TOT_NUM_FACE_VERTS];
  v->num_orig_particles = hdr[NUM_ORIG_PARTICLES];

  if (v->num_verts > 0)
    v->save_verts = new float[3 * v->num_verts];
  if (v->num_cells > 0) {
    v->num_cell_verts = new int[v->num_cells];
    v->sites = new float[3 * v->num_orig_particles];
  }
  if (v->tot_num_cell_verts > 0)
    v->cells = new int[v->tot_num_cell_verts];
  if (v->num_complete_cells > 0) {
    v->complete_cells = new int[v->num_complete_cells];
    v->areas = new float[v->num_complete_cells];
    v->vols = new float[v->num_complete_cells];
    v->num_cell_faces = new int[v->num_complete_cells];
  }
  if (v->tot_num_cell_faces > 0)
    v->num_face_verts = new int[v->tot_num_cell_faces];
  if (v->tot_num_face_verts > 0)
    v->face_verts = new int[v->tot_num_face_verts];

  struct map_block_t map[] = {

    { DIY_FLOAT,  OFST, 3, 
      offsetof(struct vblock_t, mins)                 },
    { DIY_FLOAT, ADDR, v->num_verts * 3, 
      DIY_Addr(v->save_verts)                         },
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

  DIY_Create_struct_datatype(DIY_Addr(v), 10, map, dtype);

  return v;

}
//--------------------------------------------------------------------------
//
// parse args
//
void ParseArgs(int argc, char ** argv) {

  if (argc < 13) {
    fprintf(stderr, "Usage: dense <filename> "
	    "<data size x y z> <resample grid size x y z> <swap (0 or 1)>");
    exit(0);
  }

  strcpy(infile, argv[1]);
  data_mins[0] = atoi(argv[2]);
  data_mins[1] = atoi(argv[3]);
  data_mins[2] = atoi(argv[4]);
  data_maxs[0] = atoi(argv[5]);
  data_maxs[1] = atoi(argv[6]);
  data_maxs[2] = atoi(argv[7]);
  glo_num_idx[0] = atoi(argv[8]);
  glo_num_idx[1] = atoi(argv[9]);
  glo_num_idx[2] = atoi(argv[10]);
  mass = atof(argv[11]);
  swap_bytes = atoi(argv[12]);

  // physical size of one grid step
  grid_step_size[0] = (data_maxs[0] - data_mins[0]) / (glo_num_idx[0] - 1);
  grid_step_size[1] = (data_maxs[1] - data_mins[0]) / (glo_num_idx[1] - 1);
  grid_step_size[2] = (data_maxs[2] - data_mins[0]) / (glo_num_idx[2] - 1);

}
//--------------------------------------------------------------------------
//
// grid parameters of one local block
//
// lid: block local id
// block_min_idx: global grid idx of block minimum grid point (output) (x,y,z)
// block_max_idx: global grid idx of block maximum grid point (output) (x,y,z)
// block_num_idx: number of grid points in block (output) (x,y,z)
//
void BlockGridParams(int lid, int *block_min_idx, int *block_max_idx,
		     int *block_num_idx) {

  // global grid index of block minimum grid point
  // ceilf() puts this point just inside the minimum block edge
  block_min_idx[0] = ceilf((vblocks[lid]->mins[0] - data_mins[0]) / 
				grid_step_size[0]);
  block_min_idx[1] = ceilf((vblocks[lid]->mins[1] - data_mins[1]) / 
				grid_step_size[1]);
  block_min_idx[2] = ceilf((vblocks[lid]->mins[2] - data_mins[2]) / 
				grid_step_size[2]);

  // global grid index of block maximum grid point
  // integer trunctation puts this index just indside the max block edge
  block_max_idx[0] = (vblocks[lid]->maxs[0] - data_mins[0]) / 
    grid_step_size[0];
  block_max_idx[1] = (vblocks[lid]->maxs[1] - data_mins[1]) / 
    grid_step_size[1];
  block_max_idx[2] = (vblocks[lid]->maxs[2] - data_mins[2]) / 
    grid_step_size[2];

  // eliminate duplication at the maximum block border
  if (fabs(data_mins[0] + block_max_idx[0] * grid_step_size[0] -
	   vblocks[lid]->maxs[0]) < eps && 
      fabs(vblocks[lid]->maxs[0] - data_maxs[0]) > eps)
    block_max_idx[0]--;
  if (fabs(data_mins[1] + block_max_idx[1] * grid_step_size[1] -
      vblocks[lid]->maxs[1]) < eps &&
      fabs(vblocks[lid]->maxs[1] - data_maxs[1]) > eps)
    block_max_idx[1]--;
  if (fabs(data_mins[2] + block_max_idx[2] * grid_step_size[2] -
      vblocks[lid]->maxs[2]) < eps &&
      fabs(vblocks[lid]->maxs[2] - data_maxs[2]) > eps)
    block_max_idx[2]--;

  // compute number of grid points in local block
  block_num_idx[0] = block_max_idx[0] - block_min_idx[0] + 1;
  block_num_idx[1] = block_max_idx[1] - block_min_idx[1] + 1;
  block_num_idx[2] = block_max_idx[2] - block_min_idx[2] + 1;

}
//--------------------------------------------------------------------------
//
// get cell bounds
//
// vblock: one voronoi block
// cell: current cell counter (updated by this function)
// face: current face counter (updated by this function)
// vert: current vertex counter (updated by this function)
// cell_min, cell_max: cell bounds (output)
// centroid: centroid, mean of all vertices (output)
//
void CellBounds(vblock_t *vblock, int &cell, int &face, int &vert,
		float *cell_min, float *cell_max, float *centroid) {

  centroid[0] = 0.0;
  centroid[1] = 0.0;
  centroid[2] = 0.0;
  int tot_verts = 0;

  // get cell bounds
  for (int k = 0; k < vblock->num_cell_faces[cell]; k++) { // faces

    for (int l = 0; l < vblock->num_face_verts[face]; l++) { // vertices

      int v = vblock->face_verts[vert];
	  
      if (k == 0 && l == 0 || vblock->save_verts[3 * v] < cell_min[0])
	cell_min[0] = vblock->save_verts[3 * v];
      if (k == 0 && l == 0 || vblock->save_verts[3 * v] > cell_max[0])
	cell_max[0] = vblock->save_verts[3 * v];
      centroid[0] += vblock->save_verts[3 * v];

      if (k == 0 && l == 0 || vblock->save_verts[3 * v + 1] < cell_min[1])
	cell_min[1] = vblock->save_verts[3 * v + 1];
      if (k == 0 && l == 0 || vblock->save_verts[3 * v + 1] > cell_max[1])
	cell_max[1] = vblock->save_verts[3 * v + 1];
      centroid[1] += vblock->save_verts[3 * v + 1];

      if (k == 0 && l == 0 || vblock->save_verts[3 * v + 2] < cell_min[2])
	cell_min[2] = vblock->save_verts[3 * v + 2];
      if (k == 0 && l == 0 || vblock->save_verts[3 * v + 2] > cell_max[2])
	cell_max[2] = vblock->save_verts[3 * v + 2];
      centroid[2] += vblock->save_verts[3 * v + 2];

      vert++;
      tot_verts++;

    } // vertices

    face++;

  } // faces

  centroid[0] /= tot_verts;
  centroid[1] /= tot_verts;
  centroid[2] /= tot_verts;

} 
//--------------------------------------------------------------------------
//
// grid parameters of one cell
//
// cell_grid_pts: number of grid points (not spaces) covered by cell (output) 
//  (x,y,z)
// cell_min_grid_idx: global grid index of minimum grid point in cell (output)
//  (x,y,z)
// cell_min_grid_pos: physical position of minimum grid point in cell (output)
//  (x,y,z)
// cell_mins: minimum cell vertex (x,y,z)
// cell_maxs: maximum cell vertex (x,y,z)
//
void CellGridParams(int *cell_grid_pts, int *cell_min_grid_idx,
		    float *cell_min_grid_pos, float *cell_mins, 
		    float *cell_maxs) {

  // global grid index of cell minimum grid point
  // ceilf used to put point just past minimum of cell
  cell_min_grid_idx[0] = ceilf((cell_mins[0] - data_mins[0]) / 
			       grid_step_size[0]);
  cell_min_grid_idx[1] = ceilf((cell_mins[1] - data_mins[1]) / 
			       grid_step_size[1]);
  cell_min_grid_idx[2] = ceilf((cell_mins[2] - data_mins[2]) / 
			       grid_step_size[2]);

  // global grid index of cell maximum grid point
  // integer division used to put point just short of maximum of cell
  int cell_max_grid_idx[3];
  cell_max_grid_idx[0] = (cell_maxs[0] - data_mins[0]) / grid_step_size[0];
  cell_max_grid_idx[1] = (cell_maxs[1] - data_mins[1]) / grid_step_size[1];
  cell_max_grid_idx[2] = (cell_maxs[2] - data_mins[2]) / grid_step_size[2];

  // cell minimum grid point physical position
  phys_pos(cell_min_grid_idx, cell_min_grid_pos);

  // number of grid points covered by cell bounding box
  cell_grid_pts[0] = cell_max_grid_idx[0] - cell_min_grid_idx[0] + 1;
  cell_grid_pts[1] = cell_max_grid_idx[1] - cell_min_grid_idx[1] + 1;
  cell_grid_pts[2] = cell_max_grid_idx[2] - cell_min_grid_idx[2] + 1;
  if (cell_min_grid_pos[0] > cell_maxs[0] || 
      cell_min_grid_pos[1] > cell_maxs[1] ||
      cell_min_grid_pos[2] > cell_maxs[2]) {
    cell_grid_pts[0] = 0;
    cell_grid_pts[1] = 0;
    cell_grid_pts[2] = 0;
  }

}
//--------------------------------------------------------------------------
//
// grid bounds of one cell
//
// cell_mins: minimum cell vertex (x,y,z)
// cell_maxs: maximum cell vertex (x,y,z)
// block_min_idx, block_max_idx: min/max block global grid point indices (x,y,z)
// num_cell_grid_pts: number of grid points (not spaces) covered by cell 
//  (output) (x,y,z)
// cell_min_grid_idx: global grid index of minimum grid point in cell 
//  (output) (x,y,z)
// cell_max_grid_idx: global grid index of maximum grid point in cell
//  (output) (x,y,z)
//
void CellGridBounds(float *cell_mins, float *cell_maxs, int *block_min_idx,
		    int *block_max_idx, int *num_cell_grid_pts, 
		    int *cell_min_grid_idx, int *cell_max_grid_idx) {

  // global grid index of cell minimum grid point
  // ceilf used to put point just past minimum of cell
  cell_min_grid_idx[0] = (cell_mins[0] - data_mins[0]) / grid_step_size[0];
  cell_min_grid_idx[1] = (cell_mins[1] - data_mins[1]) / grid_step_size[1];
  cell_min_grid_idx[2] = (cell_mins[2] - data_mins[2]) / grid_step_size[2];

  // global grid index of cell minimum grid point
  // integer division used to put point just short of maximum of cell
  cell_max_grid_idx[0] = (cell_maxs[0] - data_mins[0]) / grid_step_size[0];
  cell_max_grid_idx[1] = (cell_maxs[1] - data_mins[1]) / grid_step_size[1];
  cell_max_grid_idx[2] = (cell_maxs[2] - data_mins[2]) / grid_step_size[2];

  // number of grid steps covered by cell
  num_cell_grid_pts[0] = cell_max_grid_idx[0] - cell_min_grid_idx[0] + 1;
  num_cell_grid_pts[1] = cell_max_grid_idx[1] - cell_min_grid_idx[1] + 1;
  num_cell_grid_pts[2] = cell_max_grid_idx[2] - cell_min_grid_idx[2] + 1;

  // fixup grid point ranges that are too small
  float cell_min_grid_pos[3];
  float cell_max_grid_pos[3];
  phys_pos(cell_min_grid_idx, cell_min_grid_pos);
  phys_pos(cell_max_grid_idx, cell_max_grid_pos);
  while (num_cell_grid_pts[0] < 2) {
    if (cell_min_grid_pos[0] > cell_mins[0] && 
	cell_min_grid_idx[0] > block_min_idx[0])
      cell_min_grid_idx[0]--;
    else if (cell_max_grid_pos[0] < cell_maxs[0] && 
	     cell_max_grid_idx[0] < block_max_idx[0])
      cell_max_grid_idx[0]++;
    num_cell_grid_pts[0]++;
  }
  while (num_cell_grid_pts[1] < 2) {
    if (cell_min_grid_pos[1] > cell_mins[1] && 
	cell_min_grid_idx[1] > block_min_idx[1])
      cell_min_grid_idx[1]--;
    else if (cell_max_grid_pos[1] < cell_maxs[1] && 
	     cell_max_grid_idx[1] < block_max_idx[1])
      cell_max_grid_idx[1]++;
    num_cell_grid_pts[1]++;
  }
  while (num_cell_grid_pts[2] < 2) {
    if (cell_min_grid_pos[2] > cell_mins[2] && 
	cell_min_grid_idx[2] > block_min_idx[2])
      cell_min_grid_idx[2]--;
    else if (cell_max_grid_pos[2] < cell_maxs[2] && 
	     cell_max_grid_idx[2] < block_max_idx[2])
      cell_max_grid_idx[2]++;
    num_cell_grid_pts[2]++;
  }

}
//--------------------------------------------------------------------------
//
// distributes scalar value to grid points within a window size of an 
//  input point
//
// pt: input point
// scalar: scalar value at input point
// win_min, win_max: physical extents of desired window (x,y,z)
// block_min_idx, block_max_idx: min/max block global grid point indices (x,y,z)
// grid_idxs: global grid indices of grid points within window size of
//  input point (x,y,z,x,y,z,...) (output)
// grid_scalars: distributed scalars at each grid_idx (output)
//
void DistributeScalar(float *pt, float scalar, float *win_min, float *win_max,
		      int *block_min_idx, int *block_max_idx, 
		      vector <int> &grid_idxs, vector <float> &grid_scalars) {

  // global grid indices of window min and max grid points
  int min_win_idx[3];
  int max_win_idx[3];

  // find window bounds
  int num_win_idx[3];
  CellGridBounds(win_min, win_max, block_min_idx, block_max_idx, 
		 num_win_idx, min_win_idx, max_win_idx);

  float tot_weight = 0.0f; // total of weights in the window, should be 1.0

  // distribute fractional densities onto grid points in the window
  float v0 = 0.0; // volume of first box computed
  vector <float> weights; // weights accociated with grid points in the windonw
  int ijk[3]; // grid index
  for (ijk[2] = min_win_idx[2]; ijk[2] <= max_win_idx[2]; ijk[2]++) {
    for (ijk[1] = min_win_idx[1]; ijk[1] <= max_win_idx[1]; ijk[1]++) {
      for (ijk[0] = min_win_idx[0]; ijk[0] <= max_win_idx[0]; ijk[0]++) {

	grid_idxs.push_back(ijk[0]);
	grid_idxs.push_back(ijk[1]);
	grid_idxs.push_back(ijk[2]);

	// move point a little if it lies on a grid line
	float grid_pos[3]; // physical position of grid point
	float p[3] = {pt[0], pt[1], pt[2]}; // temp copy of pt
	phys_pos(ijk, grid_pos);
	if (fabs(p[0] - grid_pos[0]) < eps)
	  p[0] += 2 * eps;
	if (fabs(p[1] - grid_pos[1]) < eps) 
	  p[1] += 2 * eps;
	if (fabs(p[2] - grid_pos[2]) < eps)
	  p[2] += 2 * eps;

	// volume of box formed by input point and grid point
	float vol = fabs((grid_pos[0] - p[0]) * (grid_pos[1] - p[1]) *
			 (grid_pos[2] - p[2]));
	assert(vol > 0.0f); // sanity
	if (v0 == 0.0) // set v0 to first volume computed
	  v0 = vol;

	// debug
// 	if (vol == 0.0)
// 	  fprintf(stderr, "vol = %.3e\n grid_pos = [%.3f %.3f %.3f] "
// 		  "p = [%.3f %.3f %.3f]", vol, 
// 		  grid_pos[0], grid_pos[1], grid_pos[2],
// 		  p[0], p[1], p[2]);

	float v = v0 / vol; // volume as a factor of v0
	weights.push_back(v);
	tot_weight += v;

      }
    }
  }

  // normalize weights and deposit densities
  float tot_norm_weight = 0.0f; // total normalized weight
  for (int i = 0; i < (int)weights.size(); i++) {
    weights[i] /= tot_weight; // normalized weight
    grid_scalars.push_back(weights[i] * scalar); // scalar on the grid point
    tot_norm_weight += weights[i]; // for sanity check later, add to 1.0
  }

  // debug
//   if (num_win_idx[0] < 2 || num_win_idx[1] < 2 || num_win_idx[2] < 2)
//     fprintf(stderr, "num_win_grid_points = %d total [%d %d %d] "
// 	    "win_grid_min = [%d %d %d] "
// 	    "win_grid_max = [%d %d %d]\n", weights.size(), 
// 	    num_win_idx[0], num_win_idx[1], num_win_idx[2],
// 	    min_win_idx[0], min_win_idx[1], min_win_idx[2],
// 	    max_win_idx[0], max_win_idx[1], max_win_idx[2]);

  // sanity checks
  assert(num_win_idx[0] >= 2 && num_win_idx[1] >= 2 && num_win_idx[2] >= 2);

  // debug
//   fprintf(stderr, "tot_norm_weight = %.3f\n", tot_norm_weight);

  assert(fabs(tot_norm_weight - 1.0f) < eps); // sanity

}
//--------------------------------------------------------------------------
// //
// // DEPRECATED
// //
// //
// // distributes density to grid points within a window size of an input point
// //
// // pt: input point
// // density: density at input point
// // win_frac: window size per side as a fraction of one grid space
// //  eg. 1.0 is a window the size of +/- one grid space in x, y, z 
// //  centered on the input point (box window, not spherical)
// // grid_idxs: global grid indices of grid points within window size of
// //  input point (x,y,z,x,y,z,...) (output)
// // grid_denses: distributed densities at each grid_idx
// //
// void DistributeDensity(float *pt, float dense, float win_frac, 
// 		       vector <int> &grid_idxs, vector <float> &grid_denses) {

//   // global grid indices of window min and max grid points
//   int min_win_idx[3];
//   int max_win_idx[3];
//   float win_size[3]; // physical window size in x, y, z

//   win_size[0] = win_frac * grid_step_size[0];
//   win_size[1] = win_frac * grid_step_size[1];
//   win_size[2] = win_frac * grid_step_size[2];

//   // global grid index of window minimum grid point
//   // ceilf used to put point just past minimum of cell
//   min_win_idx[0] = ceilf((pt[0] - win_size[0]) / grid_step_size[0]);
//   min_win_idx[1] = ceilf((pt[1] - win_size[1]) / grid_step_size[1]);
//   min_win_idx[2] = ceilf((pt[2] - win_size[2]) / grid_step_size[2]);

//   // global grid index of window maximum grid point
//   // integer division used to put point just short of maximum of cell
//   max_win_idx[0] = (pt[0] + win_size[0]) / grid_step_size[0];
//   max_win_idx[1] = (pt[1] + win_size[1]) / grid_step_size[1];
//   max_win_idx[2] = (pt[2] + win_size[2]) / grid_step_size[2];

//   // total number of grid points in window
//   int num_win_idx = (max_win_idx[0] - min_win_idx[0] + 1) *
//     (max_win_idx[1] - min_win_idx[1] + 1) *
//     (max_win_idx[2] - min_win_idx[2] + 1);

//   // total volume of window
//   float min_win_pos[3];
//   float max_win_pos[3];
//   phys_pos(min_win_idx, min_win_pos);
//   phys_pos(max_win_idx, max_win_pos);
//   float win_vol = (max_win_pos[0] - min_win_pos[0]) *
//     (max_win_pos[1] - min_win_pos[1]) * (max_win_pos[2] - min_win_pos[2]);

//   float tot_weight = 0.0f; // total of weights in the window, should be 1.0

//   // distribute fractional densities onto grid points in the window
//   int ijk[3]; // grid index
//   for (ijk[2] = min_win_idx[2]; ijk[2] <= max_win_idx[2]; ijk[2]++) {
//     for (ijk[1] = min_win_idx[1]; ijk[1] <= max_win_idx[1]; ijk[1]++) {
//       for (ijk[0] = min_win_idx[0]; ijk[0] <= max_win_idx[0]; ijk[0]++) {

// 	grid_idxs.push_back(ijk[0]);
// 	grid_idxs.push_back(ijk[1]);
// 	grid_idxs.push_back(ijk[2]);

// 	// move point a little if it lies on a grid line
// 	float grid_pos[3]; // physical position of grid point
// 	phys_pos(ijk, grid_pos);
// 	if (fabs(pt[0] - grid_pos[0]) < eps) {
// 	  if (ijk[0] == min_win_idx[0])
// 	    pt[0] += eps;
// 	  else
// 	    pt[0] -= eps;
// 	}
// 	if (fabs(pt[1] - grid_pos[1]) < eps) {
// 	  if (ijk[1] == min_win_idx[1])
// 	    pt[1] += eps;
// 	  else
// 	    pt[1] -= eps;
// 	}
// 	if (fabs(pt[2] - grid_pos[2]) < eps) {
// 	  if (ijk[2] == min_win_idx[2])
// 	    pt[2] += eps;
// 	  else
// 	    pt[2] -= eps;
// 	}

// 	// volume of box formed by input point and grid point
// 	float vol = fabs((grid_pos[0] - pt[0]) * (grid_pos[1] - pt[1]) *
// 			 (grid_pos[2] - pt[2]));
// 	float weight = (win_vol - vol) / win_vol / (num_win_idx - 1);
// 	tot_weight += weight;

// 	// debug
// // 	fprintf(stderr, "weight = %.3f tot_weight = %.3f\n", weight, tot_weight);

// 	// todo: have not checked if this works for > 8 points, doubt it
// 	// for now, check if this happens and print a warning
// 	if (num_win_idx != 8)
// 	  fprintf(stderr, "Warning: number of grid points in the window is %d. "
// 		  "DistributeDensity() has only been written and tested "
// 		  "for 8 grid points\n", num_win_idx);

// 	// fractional density on the grid point
// 	grid_denses.push_back(weight * dense);

//       }
//     }
//   }

//   // debug
// //   fprintf(stderr, "total weight in the window is %.3f, should be 1.0\n", 
// // 	  tot_weight);

//   assert(fabs(tot_weight - 1.0f) < eps); // sanity

// }
//--------------------------------------------------------------------------
//
// parameters of one grid point in the bounding box of a cell
//
// grid_idx: global grid index of the point (output) (x,y,z)
// block_grid_idx: local grid index in block of the point (output) (x,y,z)
// grid_pos: physical position of the point (output) (x,y,z)
// xi, yi, zi: indices of grid point in the bounding box of the cell
// cell_min_grid_idx: global grid index of minimum grid point in cell (output)
//  (x,y,z)
// block_min_grid_idx: global grid idx of block minimum grid point (output) 
//  (x,y,z)
// cell_min_grid_pos: physical position of minimum grid point in cell (output)
//  (x,y,z)
//
// returns: whether point is in interior of cell (1 = interior, 0 = exterior)
//
int CellGridPointParams(int *grid_idx, int *block_grid_idx, float *grid_pos,
			int xi, int yi, int zi, int *cell_min_grid_idx,
			int *block_min_grid_idx, float *cell_min_grid_pos) {

  // compute global index of current grid point
  grid_idx[0] = cell_min_grid_idx[0] + xi;
  grid_idx[1] = cell_min_grid_idx[1] + yi;
  grid_idx[2] = cell_min_grid_idx[2] + zi;

  // compute local index of current grid point in this block
  block_grid_idx[0] = grid_idx[0] - block_min_grid_idx[0];
  block_grid_idx[1] = grid_idx[1] - block_min_grid_idx[1];
  block_grid_idx[2] = grid_idx[2] - block_min_grid_idx[2];

  // compute physical position of current grid point
  // will be needed to check for interior of cell (todo)
  grid_pos[0] = cell_min_grid_pos[0] + xi * grid_step_size[0];
  grid_pos[1] = cell_min_grid_pos[1] + yi * grid_step_size[1];
  grid_pos[2] = cell_min_grid_pos[2] + zi * grid_step_size[2];

  // todo: check for interior; for now just return 1
  return 1;

}
//--------------------------------------------------------------------------
//
// makes DIY datatype for sending and receiving one item
//
// dtype: pointer to the datatype
//
void ItemDtype(DIY_Datatype *dtype) {

  struct map_block_t map[] = {
    {DIY_INT,   OFST, 3, offsetof(grid_pt_t, idx)}, // global grid index
    {DIY_FLOAT, OFST, 1, offsetof(grid_pt_t, dense)}, // density
  };
  DIY_Create_struct_datatype(0, 2, map, dtype);

}
//--------------------------------------------------------------------------
//
// write density grid
//
// density: density field
// comm: MPI communicator
// nblocks: local number of blocks
// mblocks: max number of blocks in any process
//
void WriteGrid(float **density, MPI_Comm comm, int nblocks, int mblocks) {

  MPI_Status status;
  int pts_written;
  MPI_File fd; 
  int block_min_idx[3]; // global grid index of block minimum grid point
  int block_max_idx[3]; // global grid index of block maximum grid point
  int block_num_idx[3]; // number of grid points in local block
  int sizes[3]; // sizes of global array
  int subsizes[3]; // sizes of subarrays
  int starts[3]; // starting offsets of subarrays
  MPI_Datatype dtype; // subarray datatype

  // open
  int retval = MPI_File_open(comm, (char *)outfile,
			     MPI_MODE_WRONLY | MPI_MODE_CREATE,
			     MPI_INFO_NULL, &fd);
  assert(retval == MPI_SUCCESS);
  MPI_File_set_size(fd, 0); // start with an empty file every time

  // write
  for (int block = 0; block < mblocks; block++) {

    if (block < nblocks) { // non-null block

      // get local block grid parameters
      BlockGridParams(block, block_min_idx, block_max_idx, block_num_idx);

      // reversed order intentional
      sizes[0] = glo_num_idx[2];
      sizes[1] = glo_num_idx[1];
      sizes[2] = glo_num_idx[0];
      starts[0] = block_min_idx[2];
      starts[1] = block_min_idx[1];
      starts[2] = block_min_idx[0];
      subsizes[0] = block_num_idx[2];
      subsizes[1] = block_num_idx[1];
      subsizes[2] = block_num_idx[0];

      MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_C,
			       MPI_FLOAT, &dtype);
      MPI_Type_commit(&dtype);
      MPI_File_set_view(fd, 0, MPI_FLOAT, dtype, (char *)"native", 
			MPI_INFO_NULL);

      int num_pts = block_num_idx[0] * block_num_idx[1] * block_num_idx[2];

      // write block
      int errcode = MPI_File_write_all(fd, density[block], num_pts, 
					  MPI_FLOAT, &status);
      if (errcode != MPI_SUCCESS)
	handle_error(errcode, (char *)"MPI_File_write_all nonempty datatype", 
		     comm);
      MPI_Get_count(&status, MPI_FLOAT, &pts_written);
      assert(pts_written == num_pts);

      MPI_Type_free(&dtype);

    }

    else { // null block
      float unused;
      MPI_File_set_view(fd, 0, MPI_FLOAT, MPI_FLOAT, (char *)"native", 
			MPI_INFO_NULL);
      MPI_File_write_all(fd, &unused, 0, MPI_FLOAT, &status);
    }

  }

  // close
  MPI_File_close(&fd);

}
//--------------------------------------------------------------------------
//
// MPI error handler
// decodes and prints MPI error messages
//
void handle_error(int errcode, char *str, MPI_Comm comm) {

  char msg[MPI_MAX_ERROR_STRING];
  int resultlen;
  MPI_Error_string(errcode, msg, &resultlen);
  fprintf(stderr, "%s: %s\n", str, msg);
  MPI_Abort(comm, 1);

}
//-----------------------------------------------------------------------
//
// compute 1-d index in a block
// points in a block are listed in row major order
//
// block_grid_idx: 3d index in this block (x,y,z)
// block_num_idx: number of pts in each dimension in this block (x,y,z)
//
//
// returns: 1-d index
//
int index(int *block_grid_idx, int *block_num_idx) {

  return (block_grid_idx[2] * block_num_idx[1] * block_num_idx[0]
	  + block_grid_idx[1] * block_num_idx[0] + 
	  block_grid_idx[0]);

}
//-----------------------------------------------------------------------
//
// checks whether a grid point is in a block
//
// idx: target grid point (x,y,z)
// min_idx, max_idx: range of grid points in the block (x,y,z)
//
// returns: true / false
//
bool in_block(int *idx, int *min_idx, int *max_idx) {

  if (idx[0] >= min_idx[0] && idx[0] <= max_idx[0] &&
      idx[1] >= min_idx[1] && idx[1] <= max_idx[1] &&
      idx[2] >= min_idx[2] && idx[2] <= max_idx[2])
    return true;
  else
    return false;

}
//-----------------------------------------------------------------------
//
// physical position of a global grid point
//
// grid_idx: global grid index (x,y,z)
// grid_pos: physical position (x,y,z) (output)
//
void phys_pos(int *grid_idx, float *grid_pos) {

  grid_pos[0] = data_mins[0] + grid_idx[0] * grid_step_size[0];
  grid_pos[1] = data_mins[1] + grid_idx[1] * grid_step_size[1];
  grid_pos[2] = data_mins[2] + grid_idx[2] * grid_step_size[2];

}
//-----------------------------------------------------------------------
//
// compute normal of a face using Newell's method
//
// Newell's method is more robust than simply computing the cross product of
//   three points when the points are colinear or slightly nonplanar. 
//
// verts: input vertices (x,y,z,x,y,z,...)
// num_verts: number of (x,y,z) vertices
// normal: (output) normal, allocated by caller
//
// winding order is assumed to be consistent; for the tess data model this
//  produces an outward normal
//
void ComputeNormal(float *verts, int num_verts, float *normal) {

  normal[0] = 0.0;
  normal[0] = 0.0;
  normal[0] = 0.0;

  for (int i = 0; i < num_verts; i++) {
    int cur = i;
    int next = (i + 1) % num_verts;
    normal[0] += (verts[3 * cur + 1] - verts[3 * next + 1]) * 
      (verts[3 * cur + 2] + verts[3 * next + 2]);
    normal[1] += (verts[3 * cur + 2] - verts[3 * next + 2]) * 
      (verts[3 * cur] + verts[3 * next]);
    normal[2] += (verts[3 * cur] - verts[3 * next]) * 
      (verts[3 * cur + 1] + verts[3 * next + 1]);
  }

  // normalize
  float mag = sqrt(normal[0] * normal[0] + normal[1] * normal[1] +
		   normal[2] * normal[2]);
  normal[0] /= mag;
  normal[1] /= mag;
  normal[2] /= mag;

  // direction is inward, need to invert
  normal[0] *= -1.0;
  normal[1] *= -1.0;
  normal[2] *= -1.0;

}
//--------------------------------------------------------------------------
//
// compute distance from point to plane
//
// pt: point
// verts: input vertices of polygon (x,y,z,x,y,z,...)
// num_verts: number of (x,y,z) vertices
//
// returns: signed distance from point to plane
//
float PtPlaneDist(float *pt, float *verts, int num_verts) {

  float norm[3]; // normal to polygon
  float v[3]; // p - first vertex

  // normal
  ComputeNormal(verts, num_verts, norm);

  // p - first vertex
  v[0] = pt[0] - verts[0];
  v[1] = pt[1] - verts[1];
  v[2] = pt[2] - verts[2];

  // dot product between normal and v
  return(norm[0] * v[0] + norm[1] * v[1] + norm[2] * v[2]);

}
//--------------------------------------------------------------------------
//
// whether a point lies inside a cell
//
// pt: point
// vblock: one voronoi block
// cell: current cell counter (updated by this function)
// face: current face counter (updated by this function)
// vert: current vertex counter (updated by this function)
//
// returns whether point is in cell (true) or not (false)
//
bool PtInCell(float *pt, vblock_t *vblock, int &cell, int &face, int &vert) {

  int sign; // sign of distance (1 or -1)
  int old_sign = 0; // previous sign, 0 = uninitialized
  float dist; // signed distance from point to plane

  for (int k = 0; k < vblock->num_cell_faces[cell]; k++) { // faces

    // compute distance from point to face
    dist =  PtPlaneDist(pt, &(vblock->save_verts[vblock->face_verts[vert]]), 
			vblock->num_face_verts[face]);

    // check sign of distance
    sign = (dist >= 0 ? 1 : -1);
    if (old_sign == 0)
      old_sign = sign;
    if (old_sign != sign)
      return false;

    vert += vblock->num_face_verts[face]; // set vert to start of next face
    face++;

  }

  return true;

}
//--------------------------------------------------------------------------
//
// grid points covered by one cell
//
//  if the cell covers at least one grid point, then actual number of grid points
//    will be returned and the cell mass will be distributed evenly  over that 
//    nubmer of points
//  if the cell does not cover any grid points, then eight points will be returned
//    and the mass of the cell will be distributed over will be distributed onto
//    those points weighted by volume fraction
//
// block: local block number
// vblock: one voronoi block
// cell: current cell counter (updated by this function)
// face: current face counter (updated by this function)
// vert: current vertex counter (updated by this function)
// cell_mins: minimum cell vertex (x,y,z)
// cell_maxs: maximum cell vertex (x,y,z)
// grid_pts: (output) grid points covered
//
// returns: number of grid points covered by this cell
//
int CellGridPts(int block, vblock_t *vblock, int &cell, int &face, int &vert,
		float *cell_mins, float *cell_maxs, 
		vector <grid_pt_t> &grid_pts) {

  // find grid points covered by bounding box of cell
  // todo: handle the case where zero points are covered

  // global grid index of cell minimum grid point
  // ceilf used to put point just past minimum of cell
  int cell_min_grid_idx[3];
  cell_min_grid_idx[0] = ceilf((cell_mins[0] - data_mins[0]) / 
			       grid_step_size[0]);
  cell_min_grid_idx[1] = ceilf((cell_mins[1] - data_mins[1]) / 
			       grid_step_size[1]);
  cell_min_grid_idx[2] = ceilf((cell_mins[2] - data_mins[2]) / 
			       grid_step_size[2]);

  // global grid index of cell minimum grid point
  // integer division used to put point just short of maximum of cell
  int cell_max_grid_idx[3];
  cell_max_grid_idx[0] = (cell_maxs[0] - data_mins[0]) / grid_step_size[0];
  cell_max_grid_idx[1] = (cell_maxs[1] - data_mins[1]) / grid_step_size[1];
  cell_max_grid_idx[2] = (cell_maxs[2] - data_mins[2]) / grid_step_size[2];

  // cell minimum grid point physical position
  float cell_min_grid_pos[3];
  phys_pos(cell_min_grid_idx, cell_min_grid_pos);

  // number of grid points covered by cell bounding box
  int cell_grid_pts[3];
  cell_grid_pts[0] = cell_max_grid_idx[0] - cell_min_grid_idx[0] + 1;
  cell_grid_pts[1] = cell_max_grid_idx[1] - cell_min_grid_idx[1] + 1;
  cell_grid_pts[2] = cell_max_grid_idx[2] - cell_min_grid_idx[2] + 1;
  if (cell_min_grid_pos[0] > cell_maxs[0] || 
      cell_min_grid_pos[1] > cell_maxs[1] ||
      cell_min_grid_pos[2] > cell_maxs[2]) {
    cell_grid_pts[0] = 0;
    cell_grid_pts[1] = 0;
    cell_grid_pts[2] = 0;
  }

  // iterate over grid points in the bounding box of the cell
  // and determine which are interior to the cell
  for (int zi = 0; zi < cell_grid_pts[2]; zi++) { // z
    for (int yi = 0; yi < cell_grid_pts[1]; yi++) { // y
      for (int xi = 0; xi < cell_grid_pts[0]; xi++) { // x

	// compute parameters of the current grid point in this cell
	int grid_idx[3]; // global index of current grid point
	int block_grid_idx[3]; // local index of current grid point 
	// in this block
	float grid_pos[3]; // physical position of current grid point

	// compute global index of current grid point
	grid_idx[0] = cell_min_grid_idx[0] + xi;
	grid_idx[1] = cell_min_grid_idx[1] + yi;
	grid_idx[2] = cell_min_grid_idx[2] + zi;

	// compute local index of current grid point in this block
	Global2LocalIdx(grid_idx, block_grid_idx, block);

	// compute physical position of current grid point
	// will be needed to check for interior of cell (todo)
	grid_pos[0] = cell_min_grid_pos[0] + xi * grid_step_size[0];
	grid_pos[1] = cell_min_grid_pos[1] + yi * grid_step_size[1];
	grid_pos[2] = cell_min_grid_pos[2] + zi * grid_step_size[2];

	if (PtInCell(grid_pos, vblock, cell, face, vert)) {

	  grid_pt_t grid_pt;
	  grid_pt.idx[0] = grid_idx[0];
	  grid_pt.idx[1] = grid_idx[1];
	  grid_pt.idx[2] = grid_idx[2];
	  // density initialized constant for all grid pts, distributed later
	  grid_pt.dense = 1.0 / vblock->vols[cell];

	}

      }
    }
  }

  // now that we have a total number of points, iterate one more time and set the
  // weight for each
  int num_pts = (int)grid_pts.size();
  for (int i = 0; i < num_pts; i++)
    grid_pts[i].dense /= (float)num_pts;


  // todo: distribute density of cells that don't hit a grid point by
  // volume fraction (distance to grid pont)?

  return num_pts;

}
//--------------------------------------------------------------------------
//
// convert global grid index to local block grid index
//
// global_idx: global grid index (i,j,k)
// local_idx: (output) local grid index in block (i,j,k)
// block: local block number
//
void Global2LocalIdx(int *global_idx, int *local_idx, int block) {

  // local block grid parameters
  int block_min_idx[3]; // global grid index of block minimum grid point
  int block_max_idx[3]; // global grid index of block maximum grid point
  int block_num_idx[3]; // number of grid points in local block

  // get local block grid parameters
  BlockGridParams(block, block_min_idx, block_max_idx, block_num_idx);

  // compute local index of current grid point in this block
  local_idx[0] = global_idx[0] - block_min_idx[0];
  local_idx[1] = global_idx[1] - block_min_idx[1];
  local_idx[2] = global_idx[2] - block_min_idx[2];

}
//--------------------------------------------------------------------------
