#include "VirtualGrid.h"

#include "SimpleMesh.h"

// C/C++ includes
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <sstream>

namespace cosmologytools
{


VirtualGrid::VirtualGrid()
{
  this->Bounds[0] = this->Bounds[1] =
  this->Bounds[2] = this->Bounds[3] =
  this->Bounds[4] = this->Bounds[5] = 0.0;

  this->Dimensions[0] =
  this->Dimensions[1] =
  this->Dimensions[2] = 10;

  this->Spacing[0] =
  this->Spacing[1] =
  this->Spacing[2] = 0.5;
}

//-----------------------------------------------------------------------------
VirtualGrid::~VirtualGrid()
{
  for(unsigned int i=0; i < this->Buckets.size(); ++i)
    {
    this->Buckets[ i ].clear();
    }
  this->Buckets.clear();
}

//-----------------------------------------------------------------------------
void VirtualGrid::RegisterMesh(SimpleMesh &M)
{
  if( M.Empty() )
    {
    return;
    }

  // STEP 0: Compute mesh bounds
  M.GetMeshBounds( this->Bounds );

  // STEP 1: Compute virtual grid spacing
  for( int i=0; i < 3; ++i )
    {
    REAL distance     = this->Bounds[i*2+1]-this->Bounds[i*2];
    this->Spacing[i] =
      static_cast<REAL>(distance/static_cast<REAL>(this->Dimensions[i]));
    } // END for all dimensions

//  std::cout << "VIRTUAL GRID PARAMETERS:";
//  std::cout << "\nX-BOUNDS: " << this->Bounds[0] << " - " << this->Bounds[1];
//  std::cout << "\nY-BOUNDS: " << this->Bounds[2] << " - " << this->Bounds[3];
//  std::cout << "\nZ-BOUNDS: " << this->Bounds[4] << " - " << this->Bounds[5];
//  std::cout << "\nX-SPACING:" << this->Spacing[0];
//  std::cout << "\nY-SPACING: " << this->Spacing[1];
//  std::cout << "\nZ-SPACING: " << this->Spacing[2];
//  std::cout << "\nDIMS:" << this->Dimensions[0] << " ";
//  std::cout << this->Dimensions[1] << " " << this->Dimensions[2];
//  std::cout << std::endl;
//  std::cout.flush();

  // STEP 2: Allocate buckets
  INTEGER ext[6];
  this->GetExtent( ext );
  this->Buckets.resize(ExtentUtilities::ComputeNumberOfCells(ext));

  // STEP 3: Inject tets to bucket
  std::vector<INTEGER> cellIds;
  for(INTEGER cellIdx=0; cellIdx < M.GetNumberOfCells(); ++cellIdx)
    {
    if( !this->InjectCell(cellIdx, M) )
      {
      std::cerr << "WARNING: cell outside of VirtualGrid bounds!\n";
      }
    } // END for all cells
}

//-----------------------------------------------------------------------------
void VirtualGrid::RegisterMesh(SimpleMesh &M, REAL spacing[3], REAL bounds[6])
{
  if( M.Empty() )
    {
    return;
    }

  // STEP 0: set bounds
  for(int i=0; i < 6; ++i)
    {this->Bounds[i] = bounds[i];}

  // STEP 1:Set virtual grid spacing
  for( int i=0; i < 3; ++i )
    {
    this->Spacing[i] = spacing[i];
    } // END for all dimensions

 // std::cout << "VIRTUAL GRID PARAMETERS:";
 // std::cout << "\nX-BOUNDS: " << this->Bounds[0] << " - " << this->Bounds[1];
 // std::cout << "\nY-BOUNDS: " << this->Bounds[2] << " - " << this->Bounds[3];
 // std::cout << "\nZ-BOUNDS: " << this->Bounds[4] << " - " << this->Bounds[5];
 // std::cout << "\nX-SPACING:" << this->Spacing[0];
 // std::cout << "\nY-SPACING: " << this->Spacing[1];
 // std::cout << "\nZ-SPACING: " << this->Spacing[2];
 // std::cout << "\nDIMS:" << this->Dimensions[0] << " ";
 // std::cout << this->Dimensions[1] << " " << this->Dimensions[2];
 // std::cout << std::endl;
 // std::cout.flush();

  // STEP 2: Allocate buckets
  INTEGER ext[6];
  this->GetExtent( ext );
  this->Buckets.resize(ExtentUtilities::ComputeNumberOfCells(ext));

  // STEP 3: Inject tets to bucket
  std::vector<INTEGER> cellIds;
  for(INTEGER cellIdx=0; cellIdx < M.GetNumberOfCells(); ++cellIdx)
    {
    if( !this->InjectCell(cellIdx, M) )
      {
      std::cerr << "WARNING: cell outside of VirtualGrid bounds!\n";
      }
    } // END for all cells
}

//-----------------------------------------------------------------------------
#define PROCESSNODE(V,min,max) {                          \
  INTEGER nodeijk[3];                                      \
  if( this->FindBucket(V,nodeijk)) {                      \
    for(int dim=0; dim < 3; ++dim) {                      \
      if( nodeijk[dim] < min[dim])                        \
        min[dim] = nodeijk[dim];                          \
      else if( nodeijk[dim] > max[dim])                  \
        max[dim] = nodeijk[dim];                          \
    } /* END for */                                       \
  }  /* END if */                                         \
  else {                                                 \
    std::cerr << "Cell node outside of virtual grid!\n"; \
    return false;                                       \
  }                                                      \
}

bool VirtualGrid::InjectCell(INTEGER cellIdx, SimpleMesh &M)
{
  // Nodes of the cell
  REAL V0[3]; REAL V1[3]; REAL V2[3]; REAL V3[3];

  // ijkmin and ijkmax store the cell extent on the virtual grid that the
  // cell with the cellIdx w.r.t. M covers.
  INTEGER ijkmin[3];
  ijkmin[0] = ijkmin[1] = ijkmin[2] = std::numeric_limits<INTEGER>::max();
  INTEGER ijkmax[3];
  ijkmax[0] = ijkmax[1] = ijkmax[2] = std::numeric_limits<INTEGER>::min();

  // STEP 0: Find mesh cell extent, i.e., ijkmin,ijkmax
  switch( M.Stride )
    {
    case 3:
      M.GetTriangleNodes(cellIdx,V0,V1,V2);
      PROCESSNODE(V0,ijkmin,ijkmax);
      PROCESSNODE(V1,ijkmin,ijkmax);
      PROCESSNODE(V2,ijkmin,ijkmax);
      break;
    case 4:
      M.GetTetNodes(cellIdx,V0,V1,V2,V3);
      PROCESSNODE(V0,ijkmin,ijkmax);
      PROCESSNODE(V1,ijkmin,ijkmax);
      PROCESSNODE(V2,ijkmin,ijkmax);
      PROCESSNODE(V3,ijkmin,ijkmax);
      break;
    default:
      std::cerr << "ERROR: Only tet and triangle meshes are supported!\n";
      return false;
    }

  // STEP 1: Inject cell
  INTEGER ext[6];
  this->GetCellExtent( ext );

  INTEGER ijk[3];
  for(INTEGER i=ijkmin[0]; i <= ijkmax[0]; ++i)
    {
    for(INTEGER j=ijkmin[1]; j <= ijkmax[1]; ++j)
      {
      for(INTEGER k=ijkmin[2]; k <= ijkmax[2]; ++k)
        {
        ijk[0]=i; ijk[1]=j; ijk[2]=k;
        if( this->WithinCellExtent(ijk))
          {
          INTEGER idx = ExtentUtilities::GetLinearIndex(i,j,k,ext);
          assert("pre: bucket index out-of-bounds" &&
              (idx >= 0) && (idx < this->Buckets.size() ) );
          this->Buckets[ idx ].insert( cellIdx );
          } // END if within cell extent
        } // END for all k
      } // END for all j
    } // END for all i
  return true;
}

//-----------------------------------------------------------------------------
std::string VirtualGrid::ToLegacyVtkString()
{
  std::ostringstream oss;
  oss.str("");
  oss << "# vtk DataFile Version 3.0\n";
  oss << "vtk output\n";
  oss << "ASCII\n";
  oss << "DATASET STRUCTURED_POINTS\n";
  oss << "DIMENSIONS " << this->Dimensions[0] << " ";
  oss << this->Dimensions[1] << " " << this->Dimensions[2] << std::endl;
  oss << "SPACING "<< this->Spacing[0] << " ";
  oss << this->Spacing[1] << " " << this->Spacing[2] << std::endl;
  oss << "ORIGIN " << this->Bounds[0] << " " << this->Bounds[2];
  oss << " " << this->Bounds[4] << std::endl;
  return( oss.str() );
}

//-----------------------------------------------------------------------------
void VirtualGrid::GetCandidateCellsForPoint(
        REAL pnt[3], std::vector<INTEGER> &cells)
{
  cells.resize(0);

  // STEP 0: Find point bucket
  INTEGER ijk[3];
  if( !this->FindBucket(pnt,ijk) )
    {
    // query point is outside the virtual grid
    return;
    }

  // STEP 1: Get the linear index
  INTEGER ext[6];
  this->GetCellExtent(ext);
  INTEGER idx = ExtentUtilities::GetLinearIndex(ijk,ext);
  assert("pre: bucket index is out-of-bounds!" &&
       (idx >= 0) && (idx < static_cast<INTEGER>(this->Buckets.size())));

  // STEP 2: Get the cells in the bucket
  std::set<INTEGER>::iterator iter = this->Buckets[idx].begin();
  for(; iter != this->Buckets[idx].end(); ++iter )
    {
    cells.push_back( *iter );
    } // END for all cells in the bucket
}

//-----------------------------------------------------------------------------
bool VirtualGrid::FindBucket(REAL pnt[3], INTEGER ijk[3])
{
  // STEP 0: Get distance of point to origin
  INTEGER d[3];
  for( int i=0; i < 3; ++i)
    {
    d[i] = pnt[i]-this->Bounds[i*2];
    if( d[i] < 0.0 )
      {
      // pnt is outside virtual grid
      return false;
      }
    ijk[i] = static_cast<INTEGER>(std::floor((d[i]/this->Spacing[i])));
    } // END for all directions

  if( !this->WithinCellExtent(ijk) )
    {
    return false;
    }
  return true;
}

} /* namespace cosmologytools */
