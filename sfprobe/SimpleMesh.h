/**
 * @class SimpleMesh
 * @brief Data-structure to represent an unstructured mesh
 */
#ifndef SIMPLEMESH_H_
#define SIMPLEMESH_H_

#include "CosmologyToolsMacros.h"

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
   * @brief
   * @param pntIdx
   * @param pnt
   */
  void GetNode(INTEGER pntIdx, REAL pnt[3]);

  /**
   * @brief Returns the IDs of the cell corresponding to the given cell index.
   * @param cellIdx the index of the cell
   * @param cellIds the IDs of the cell
   */
  void GetCell( INTEGER cellIdx, std::vector<INTEGER> &cellIds);

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
