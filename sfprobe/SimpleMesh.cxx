#include "SimpleMesh.h"

#include <cassert>

namespace cosmologytools {

SimpleMesh::SimpleMesh()
{
   this->Stride=0;
}

//------------------------------------------------------------------------------
SimpleMesh::~SimpleMesh()
{
  this->Nodes.clear();
  this->Connectivity.clear();
}

//------------------------------------------------------------------------------
void SimpleMesh::GetCell(INTEGER cellIdx, std::vector<INTEGER> &cellIds)
{
  assert("pre: cellIdx is out-of-bounds!" &&
         (cellIdx >= 0) && (cellIdx < this->GetNumberOfCells()));

  cellIds.clear();
  cellIds.resize( this->Stride );
  for( INTEGER i=0; i < this->Stride; ++i )
    {
    cellIds[ i ] = this->Connectivity[cellIdx*this->Stride+i];
    }
}

//------------------------------------------------------------------------------
void SimpleMesh::GetNode(INTEGER pntIdx, REAL pnt[3])
{
  assert("pre: pntIdx is out-of-bounds!" &&
         (pntIdx >= 0) && (pntIdx < this->GetNumberOfNodes()));

  pnt[0] = this->Nodes[ pntIdx*3  ];
  pnt[1] = this->Nodes[ pntIdx*3+1];
  pnt[2] = this->Nodes[ pntIdx*3+2];
}

} /* namespace cosmologytools */
