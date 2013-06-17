/**
 * @class SimpleMesh
 * @brief Data-structure to represent an unstructured mesh
 */
#ifndef SIMPLEMESH_H_
#define SIMPLEMESH_H_

#include "CosmoToolsMacros.h"

// C/C++ includes
#include <vector>

namespace cosmologytools {

class SimpleMesh
{
public:
  SimpleMesh();
  virtual ~SimpleMesh();

  INTEGER Stride;
  std::vector<REAL> Nodes;
  std::vector<INTEGER> Connectivity;

  /**
   * @brief Returns the number of nodes in the mesh
   * @return N total number of nodes
   */
  INTEGER GetNumberOfNodes()
    {
    return( static_cast<INTEGER>(this->Nodes.size()/3) );
    }

  /**
   * @brief Returns the number of cells in the mesh
   * @return N total number of cells
   */
  INTEGER GetNumberOfCells()
    {
    return( static_cast<INTEGER>(this->Connectivity.size()/this->Stride) );
    }

  /**
   * @brief Checks if this SimpleMesh instance is an empty mesh.
   * @return Status true iff the mesh is empty, else, false.
   * @note A mesh is considered empty iff this->GetNumberOfNodes()==0 ||
   * this->GetNumberOfCells()==0.
   */
  bool Empty()
    {
    if( (this->Nodes.size()==0) || (this->Connectivity.size()==0) )
      return true;
    return false;
    }

  /**
   * @brief Computes the bounds of the mesh
   * @param bounds array storing [xmin,xmax,ymin,ymax,zmin,zmax]
   */
  void GetMeshBounds(REAL bounds[6]);

  /**
   * @brief Returns the cartesian coordinates of the node at the given index.
   * @param pntIdx the index of the node in query (in)
   * @param pnt the cartesian nodes of the point (out)
   */
  void GetNode(INTEGER pntIdx, REAL pnt[3]);

  /**
   * @brief Returns the IDs of the cell corresponding to the given cell index.
   * @param cellIdx the index of the cell
   * @param cellIds the IDs of the cell
   */
  void GetCell(INTEGER cellIdx, std::vector<INTEGER> &cellIds);

  /**
   * @brief A convenience method to return the coordinates of the triangle nodes
   * @param cellIdx the index of the triangle cell (in)
   * @param V0 xyz coordinates of node 0 (out)
   * @param V1 xyz coordinates of node 1 (out)
   * @param V2 xyz coordinates of node 2 (out)
   * @pre (cellIdx >= 0) && (cellIdx < this->GetNumberOfCells)
   * @pre (this->Stride == 3) i.e., the mesh must be triangular
   */
  void GetTriangleNodes(
      INTEGER cellIdx, REAL V0[3], REAL V1[3], REAL V2[3]);

  /**
   * @brief A convenience method to return the coordinates of the tet nodes
   * @param cellIdx the index of the tetrahedron cell (in)
   * @param V0 xyz coordinates of node 0 (out)
   * @param V1 xyz coordinates of node 1 (out)
   * @param V2 xyz coordinates of node 2 (out)
   * @param V3 xyz coordinates of node 3 (out)
   * @pre (cellIdx >= 0) && (cellIdx < this->GetNumberOfCells)
   * @pre (this->Stride == 4) i.e., the mesh must be a tetrahedral mesh
   */
  void GetTetNodes(
     INTEGER cellIdx, REAL V0[3], REAL V1[3], REAL V2[3], REAL V3[3]);

  /**
   * @brief Clears all data of this mesh.
   * @post this->Empty()==true.
   */
  void Clear()
    {
    this->Nodes.clear();
    this->Connectivity.clear();
    }

private:
  DISABLE_COPY_AND_ASSIGNMENT(SimpleMesh);
};

} /* namespace cosmologytools */
#endif /* SIMPLEMESH_H_ */
