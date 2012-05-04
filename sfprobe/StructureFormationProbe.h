/**
 * @class StructureFormationProbe
 * @brief A class the implements a tesselation-based approach to
 */
#ifndef STRUCTUREFORMATIONPROBE_H_
#define STRUCTUREFORMATIONPROBE_H_

#include "CosmologyToolsMacros.h"

namespace cosmologytools {

class StructureFormationProbe
{
public:
  StructureFormationProbe();
  virtual ~StructureFormationProbe();

  // In-line functions
  SetVector6Macro(WholeExtent,int);
  GetNSetMacro(Particles,REAL*);
  GetNSetMacro(NumParticles,int);
  GetNSetMacro(Stride,int);

  /**
   * @brief Sets the virtual grid extent for this instance.
   * @param extent the node extent for this StructureFormationProbe instance.
   * @note extent is a flat array that stores[imin imax jmin jmax kmin kmax]
   * @pre extent is expected to be a 3-D extent
   */
  void SetGridExtent(int extent[6]);

  /**
   * @brief Tesselates the grid extent.
   * @post Each voxel is split to 5 tetrahedra.
   */
  virtual void Tesselate();

  /**
   * @brief Writes the tesselation in a VTK file.
   * @param fileName the name of the file
   * @note used for debugging.
   */
  void WriteTesselation( char* fileName);

protected:
  int Extent[6];       // The grid node extent (supplied by the user)
  int WholeExtent[6];  // The entire grid extent (supplied by the user)
  REAL *Particles;     // A flat particles array (supplied by the user)
  int NumParticles;    // The total number of particles (supplied by the user)
  int Stride;          // The stide used to index the particles array


  int NumTets;        // The total number of tets
  int  *Connectivity; // tetrahedral connectivity computed internally.

  /**
   * @brief Initializes data-structures
   */
  void Initialize();

  /**
   * @brief Tesselates the extent that is assigned to this process
   */
  void TesselateGridExtent();

  /**
   * @brief Given structured coordinates in IJK space return the linear index,
   * w.r.t. the user-supplied global grid extent (node-extent)
   * @param i the i-coordinate of the node in IJK space
   * @param j the j-coordinate of the node in IJK space
   * @param k the k-coordinate of the node in IJK space
   * @return idx the global linear index of the node
   */
  int GetGlobalLinearIndex(const int i, const int j, const int k);

  /**
   * @brief Given structured coordinates in IJK space return the linear index.
   * @param i the i-coordinate of the node in IJK space
   * @param j the j-coordinate of the node in IJK space
   * @param k the k-coordinate of the node in IJK space
   * @return idx the linear index of the node
   * @note This linear index is local w.r.t. the user-supplied grid extent.
   */
  int GetLinearIndex(const int i, const int j, const int k);

  /**
   * @brief Given a (local) linear index, returns the global IJK coordinates.
   * @param idx the linear index.
   * @param i the i-coordinate of the node in IJK space (out)
   * @param j the j-coordinate of the node in IJK space (out)
   * @param k the k-coordinate of the node in IJK space (out)
   * @post (i,j,k) must be within the user-supplied extent.
   */
  void GetStructuredCoordinates(
      const int idx,int &i, int &j, int &k);

  /**
   * @brief A convenience method that checks if the extent is a 3-D extent.
   * @param extent the grid extent to check whether or not it is 3-D.
   * @return status true iff the extent corresponds to a 3-D extent, else false.
   */
  bool Is3DExtent( int extent[6] );

  /**
   * @brief Computes the number of nodes given the grid node-extent.
   * @param extent the node extent of the grid
   * @return N the total number of nodes in the grid.
   * @pre this->Is3DExtent( extent ) == true.
   * @post N > 0
   */
  int GetNumberOfNodes( int extent[6] );

  /**
   * @brief Computes the number of cells given the grid node-extent
   * @param extent the node extent of the grid
   * @return N the total number of cells in the grid.
   * @pre this->Is3DExtent( extent ) == true.
   * @post N > 0
   */
  int GetNumberOfCells( int extent[6] );

private:
  DISABLE_COPY_AND_ASSIGNMENT(StructureFormationProbe);
};

} /* namespace hacctools */
#endif /* STRUCTUREFORMATIONPROBE_H_ */
