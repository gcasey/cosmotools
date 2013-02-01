//------------------------------------------------------------------------------
//
// serial io class
//
// Tom Peterka
// Argonne National Laboratory
// 9700 S. Cass Ave.
// Argonne, IL 60439
// tpeterka@mcs.anl.gov
//
// All rights reserved. May not be used, modified, or copied
// without permission
//
//--------------------------------------------------------------------------

#if 0 // not using compression for now
#include "zlib.h"
#endif

#include "ser_io.hpp"
#include "swap.hpp"

//----------------------------------------------------------------------------
//
// reads all blocks in a file
//
// filename: input filename
// blocks: pointers to blocks (output)
// compress: whether file is compressed (optional, default = not compressed)
//
// side effects: allocates space for the new blocks
//
// returns: total number of blocks
//
int SER_IO::ReadAllBlocks(const char *filename, vblock_t** &blocks,
			  bool compress) {

  compress = compress; // quiet compiler warning

  FILE *fd;
  int64_t *ftr; // footer
  int tb; // total number of blocks

  fd = fopen(filename, "r");
  assert(fd != NULL);

  ReadFooter(fd, ftr, tb);

  blocks = new vblock_t*[tb];
  for (int i = 0; i < tb; i++) {

#if 0 // not using compression for now

    int comp_size; // size of compressed block in bytes
    int decomp_size; // size of decompressed block in bytes
    int count;

    if (compress) {

      if (i < tb - 1)
	comp_size = ftr[i + 1] - ftr[i];
      else {
	fseek(fd, 0, SEEK_END);
	comp_size = ftell(fd) - ftr[i];
      }

      unsigned char *comp_buf = new unsigned char[comp_size];
      vector<unsigned char> decomp_buf;
      fseek(fd, ftr[i], SEEK_SET);
      count = fread(comp_buf, 1, comp_size, fd);
      assert(count == comp_size);
      DecompressBlock(comp_buf, comp_size, &decomp_buf, &decomp_size);
      CopyBlock(&decomp_buf[0], blocks[i]);
      delete[] comp_buf;

    }

    else

#endif 

      ReadBlock(fd, blocks[i], ftr[i]);

  }

  fclose(fd);

  delete[] ftr;
  return tb;

}
//----------------------------------------------------------------------------
//
// reads the file footer
// footer in file is always ordered by global block id
// output footer is in the same order
//
// fd: open file
// ftr: footer data (output)
// tb: total number of blocks in the file (output)
//
// side effects: allocates ftr
//
void SER_IO::ReadFooter(FILE*& fd, int64_t*& ftr, int& tb) {

  int ofst;
  int count;
  int64_t temp;

  ofst = sizeof(int64_t);
  fseek(fd, -ofst, SEEK_END);
  count = fread(&temp, sizeof(int64_t), 1, fd); // total number of blocks
  assert(count == 1); // total number of blocks

  if (swap_bytes)
    Swap((char *)&temp, 1, sizeof(int64_t));
  tb = temp;

  if (tb > 0) {

    ftr = new int64_t[tb];
    ofst = (tb + 1) * sizeof(int64_t);
    fseek(fd, -ofst, SEEK_END);
    count = fread(ftr, sizeof(int64_t), tb, fd);
    assert(count == tb);

    if (swap_bytes)
      Swap((char *)ftr, tb, sizeof(int64_t));

  }

}
//----------------------------------------------------------------------------
//
// reads the header for one block from a file
//
// fd: open file
// hdr: allocated header data
// ofst: location in file of the header (bytes)
//
void SER_IO::ReadHeader(FILE *fd, int *hdr, int64_t ofst) {

  int count;

  fseek(fd, ofst, SEEK_SET);
  count = fread(hdr, sizeof(int), DIY_MAX_HDR_ELEMENTS, fd);
  assert(count == DIY_MAX_HDR_ELEMENTS);

  if (swap_bytes)
    Swap((char *)hdr, DIY_MAX_HDR_ELEMENTS, sizeof(int));

}
//----------------------------------------------------------------------------
//
// Copies the header for one block from a buffer in memory
//
// in_buf: input buffer location
// hdr: allocated header data
//
// returns: number of bytes copies
//
int SER_IO::CopyHeader(unsigned char *in_buf, int *hdr) {

  memcpy(hdr, in_buf, DIY_MAX_HDR_ELEMENTS * sizeof(int));

  if (swap_bytes)
    Swap((char *)hdr, DIY_MAX_HDR_ELEMENTS, sizeof(int));

  return(DIY_MAX_HDR_ELEMENTS * sizeof(int));

}
//----------------------------------------------------------------------------
//
// reads one block from a file
//
// fd: open file
// v: pointer to output block
// ofst: file file pointer to start of header for this block
//
// side-effects: allocates block
//
void SER_IO::ReadBlock(FILE *fd, vblock_t* &v, int64_t ofst) {

  // get header info
  int hdr[DIY_MAX_HDR_ELEMENTS];
  ReadHeader(fd, hdr, ofst);

  // create block
  v = new vblock_t;
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

  fread(v->mins, sizeof(float), 3, fd);
  fread(v->save_verts, sizeof(float), 3 * v->num_verts, fd);
  fread(v->sites, sizeof(float), 3 * v->num_orig_particles, fd);
  fread(v->complete_cells, sizeof(int), v->num_complete_cells, fd);
  fread(v->areas, sizeof(float), v->num_complete_cells, fd);
  fread(v->vols, sizeof(float), v->num_complete_cells, fd);
  fread(v->num_cell_faces, sizeof(int), v->num_complete_cells, fd);
  fread(v->num_face_verts, sizeof(int), v->tot_num_cell_faces, fd);
  fread(v->face_verts, sizeof(int), v->tot_num_face_verts, fd);
  fread(v->maxs, sizeof(float), 3, fd);

  if (swap_bytes) {
    Swap((char *)v->mins, 3, sizeof(float));
    Swap((char *)v->save_verts, 3 * v->num_verts, sizeof(float));
    Swap((char *)v->sites, 3 * v->num_orig_particles, sizeof(float));
    Swap((char *)v->complete_cells, v->num_complete_cells, sizeof(int));
    Swap((char *)v->areas, v->num_complete_cells, sizeof(float));
    Swap((char *)v->vols, v->num_complete_cells, sizeof(float));
    Swap((char *)v->num_cell_faces, v->num_complete_cells, sizeof(int));
    Swap((char *)v->num_face_verts, v->tot_num_cell_faces, sizeof(int));
    Swap((char *)v->face_verts, v->tot_num_face_verts, sizeof(int));
    Swap((char *)v->maxs, 3, sizeof(float));
  }

}
//----------------------------------------------------------------------------
//
// copies one block from a buffer in memory
//
// in_buf: input buffer location
// m: pointer to output block
// ofst: file file pointer to start of header for this block
//
// side-effects: allocates block
//
void SER_IO::CopyBlock(unsigned char *in_buf, vblock_t* &v) {

  int ofst = 0; // offset in buffer for next section of data to read

  // get header info
  int hdr[DIY_MAX_HDR_ELEMENTS];
  ofst += CopyHeader(in_buf, hdr);

  // create block
  v = new vblock_t;
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

    memcpy(&v->mins, in_buf + ofst, 3 * sizeof(float));
    ofst += (3 * sizeof(float));

    memcpy(v->save_verts, in_buf + ofst, 3 * v->num_verts * sizeof(float));
    ofst += (3 * v->num_verts * sizeof(float));

    memcpy(v->sites, in_buf + ofst, 3 * v->num_orig_particles * sizeof(float));
    ofst += (3 * v->num_orig_particles * sizeof(float));

    memcpy(v->complete_cells, in_buf + ofst, 
	   v->num_complete_cells * sizeof(int));
    ofst += (v->num_complete_cells * sizeof(int));

    memcpy(v->areas, in_buf + ofst, v->num_complete_cells * sizeof(float));
    ofst += (v->num_complete_cells * sizeof(float));

    memcpy(v->vols, in_buf + ofst, v->num_complete_cells * sizeof(float));
    ofst += (v->num_complete_cells * sizeof(float));

    memcpy(v->num_cell_faces, in_buf + ofst, 
	   v->num_complete_cells * sizeof(int));
    ofst += (v->num_complete_cells * sizeof(int));

    memcpy(v->num_face_verts, in_buf + ofst, 
	   v->tot_num_cell_faces * sizeof(int));
    ofst += (v->tot_num_cell_faces * sizeof(int));

    memcpy(v->face_verts, in_buf + ofst, v->tot_num_face_verts * sizeof(int));
    ofst += (v->tot_num_face_verts * sizeof(int));

    memcpy(&v->maxs, in_buf + ofst, 3 * sizeof(float));
    ofst += (3 * sizeof(float));

  if (swap_bytes) {
    Swap((char *)&v->mins, 3, sizeof(float));
    Swap((char *)v->save_verts, 3 * v->num_verts, sizeof(float));
    Swap((char *)v->sites, 3 * v->num_orig_particles, sizeof(float));
    Swap((char *)v->complete_cells, v->num_complete_cells, sizeof(int));
    Swap((char *)v->areas, v->num_complete_cells, sizeof(float));
    Swap((char *)v->vols, v->num_complete_cells, sizeof(float));
    Swap((char *)v->num_cell_faces, v->num_complete_cells, sizeof(int));
    Swap((char *)v->num_face_verts, v->tot_num_cell_faces, sizeof(int));
    Swap((char *)v->face_verts, v->tot_num_face_verts, sizeof(int));
    Swap((char *)&v->maxs, 3, sizeof(float));
  }

}
//----------------------------------------------------------------------------

#if 0 // not using compression for now

//
// block decompression
//
// in_buf: input block buffer (MPI_BYTE datatype)
// in_size: input size in bytes
// decomp_buf: decompressed buffer
// decomp_size: decompressed size in bytes (output)
//
// side effects: grows decomp_buf, user's responsibility to free it
//
void SER_IO::DecompressBlock(unsigned  char* in_buf, int in_size, 
			    vector<unsigned char> *decomp_buf, 
			    int *decomp_size) {

  unsigned int out_size; // size of output buffer
  z_stream strm; // structure used to pass info t/from zlib
  unsigned char temp_buf[CHUNK]; // temporary
  int num_out; // number of compressed bytes in one deflation
  int ret; // return value

  // allocate inflate state
  strm.zalloc = Z_NULL;
  strm.zfree = Z_NULL;
  strm.opaque = Z_NULL;
  strm.avail_in = 0;
  strm.next_in = Z_NULL;
  ret = inflateInit(&strm);
  assert(ret == Z_OK);

  *decomp_size = 0;
  int num_in_chunks = ceil((float)in_size / CHUNK); // number of input chunks

  // read the character buffer in chunks
  for (int i = 0; i < num_in_chunks; i++) { 
    strm.next_in = in_buf + i * CHUNK;
    if (i < num_in_chunks - 1) 
      strm.avail_in = CHUNK;
    else
      strm.avail_in = (in_size % CHUNK ? in_size - i * CHUNK : CHUNK);
    do { // decompress each chunk
      strm.avail_out = CHUNK;
      strm.next_out = temp_buf;
      ret = inflate(&strm, Z_NO_FLUSH);
      assert(ret != Z_STREAM_ERROR);
      num_out = CHUNK - strm.avail_out;
      *decomp_size += num_out;
      decomp_buf->insert(decomp_buf->end(), temp_buf, temp_buf + num_out);
    } while (strm.avail_out == 0);
    if (ret == Z_STREAM_END) // in case decompression finishes early
      break;
  }
  inflateEnd(&strm);

//   fprintf(stderr, "decompression: from %d bytes to %d bytes\n", 
// 	  in_size, *decomp_size);

}
//-----------------------------------------------------------------------

#endif
