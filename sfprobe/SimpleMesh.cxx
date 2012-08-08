#include "SimpleMesh.h"

#include <limits>
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

//------------------------------------------------------------------------------
void SimpleMesh::GetMeshBounds(REAL bounds[6])
{

 // Initialize bounds
 bounds[0] = std::numeric_limits<REAL>::max(); // xmin
 bounds[1] = std::numeric_limits<REAL>::min(); // xmax
 bounds[2] = std::numeric_limits<REAL>::max(); // ymin
 bounds[3] = std::numeric_limits<REAL>::min(); // ymax
 bounds[4] = std::numeric_limits<REAL>::max(); // zmin
 bounds[5] = std::numeric_limits<REAL>::min(); // zmax

 REAL pnt[3];
 for(INTEGER node=0; node < this->GetNumberOfNodes(); ++node)
   {
   this->GetNode(node,pnt);
   for( int i=0; i < 3; ++i )
     {
     if(pnt[i] < bounds[i*2])
       {
       bounds[i*2] = pnt[i];
       }
     else if(pnt[i] > bounds[i*2+1])
       {
       bounds[i*2+1] = pnt[i];
       }
     } // END for all dimensions
   } // END for all nodes
}


} /* namespace cosmologytools */
