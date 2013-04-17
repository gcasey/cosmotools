#ifndef HACC_PARTICLES_IO_H_
#define HACC_PARTICLES_IO_H_

#include <stdint.h> // For int64_t
#include <mpi.h> // For MPI_Comm

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initializes hacc io data-structures
 * @param comm MPI communicator to use (in)
 */
void hacc_io_init(MPI_Comm *comm);

/**
 * @brief Initializes hacc io data-structures with Fortran communicator
 * @param comm MPI communicator to use (in)
 * @note Fortran codes should call this method instead of init
 * @see hacc_io_init
 */
void hacc_io_finit(MPI_Comm *comm);

/**
 * @brief Sets the file handle
 * @param file the path and filename of the file to use.
 */
void hacc_io_set_file(char *file);

/**
 * @brief Returns the number of elements to read by this rank
 * @param N number of elements
 * @note Call this method to find out how big the arrays must be when
 * you call hacc_io_read_particles
 */
void hacc_io_get_num_elements(int *N);

/**
 * @brief Reads in the particles and associated data.
 * @param x the x-coordinates of all the particles of this rank.
 * @param y the y-coordinates of all the particles of this rank.
 * @param z the z-coordinates of all the particles of this rank.
 * @param vx the x velocity component of all the particles of this rank.
 * @param vy the y velocity component of all the particles of this rank.
 * @param vz the z velocity component of all the particles of this rank.
 * @param phi potential? of all the particles of this rank.
 * @param idx the global ID of all the particles of this rank.
 * @pre arrays must be pre-allocated
 */
void hacc_io_read_particles(
      float *x, float *y, float *z,
      float *vx, float *vy, float *vz,
      float *phi,
      int64_t *idx);

/**
 * @brief De-allocates all dynamically allocated objects.
 */
void hacc_io_finalize();

#ifdef __cplusplus
}
#endif

#endif /* HACC_PARTICLES_IO_H_ */
