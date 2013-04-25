//=============================================================================
//
//  Copyright (c) Kitware, Inc.
//  All rights reserved.
//  See LICENSE.txt for details.
//
//  This software is distributed WITHOUT ANY WARRANTY; without even
//  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
//  PURPOSE.  See the above copyright notice for more information.
//
//  Copyright 2012 Sandia Corporation.
//  Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
//  the U.S. Government retains certain rights in this software.
//
//=============================================================================
#ifndef __TetStructure_h
#define __TetStructure_h

#include <dax/CellTag.h>
#include <dax/cont/ArrayHandle.h>
#include <dax/cont/UnstructuredGrid.h>
#include <dax/exec/CellField.h>
#include <dax/exec/CellVertices.h>
#include <dax/exec/ExecutionObjectBase.h>
#include <dax/math/Compare.h>
#include <dax/math/Sign.h>

#include "UtilityFunctions.h"

namespace dax {
namespace exec {

class TetStructure : public dax::exec::ExecutionObjectBase
{
  //we currently don't support custom ArrayHandle for any of the portals
  //those could become template parameters of the TetStructure class
  //so that we could be more general
  typedef dax::exec::internal::TopologyUnstructured<
            dax::CellTagTetrahedron,
            dax::cont::ArrayHandle<dax::Id>::PortalConstExecution >
        TopologyPortal;
typedef dax::cont::ArrayHandle<dax::Vector3>::PortalConstExecution
        CoordinatesPortal;
  typedef dax::cont::ArrayHandle<dax::Scalar>::PortalConstExecution
        VolumePortal;

  TopologyPortal Topo;
  CoordinatesPortal Coords;
  VolumePortal Volumes;
public:

  typedef  dax::exec::CellVertices< dax::CellTagTetrahedron > CellVerts;
  typedef  dax::exec::CellField< dax::Vector3, dax::CellTagTetrahedron > CellCoords;

  DAX_CONT_EXPORT TetStructure(
     dax::cont::UnstructuredGrid<dax::CellTagTetrahedron> grid,
     dax::cont::ArrayHandle<dax::Scalar> volumes):
  Topo(grid.PrepareForInput()),
  Coords(grid.GetPointCoordinates().PrepareForInput()),
  Volumes(volumes.PrepareForInput())

  {
    //by constructing a TetStructure you transfer all the tet information
    //from the host to the device machine
  }

  DAX_EXEC_EXPORT dax::Id GetNumberOfCells( ) const
    {
    return this->Topo.GetNumberOfCells();
    }

  //Returns 1 if the point is inside the given cell
  //will increment volume with the volume of the cell if it was valid
  DAX_EXEC_EXPORT dax::Id ProbePoint(dax::Id index, dax::Vector3 point,
                                     dax::Scalar& totalVolume) const
    {
    const CellVerts verts = this->Topo.GetCellConnections( index );
    CellCoords coords;
    for (int i = 0; i < 4; i++)
      { coords[i] = this->Coords.Get(verts[i]); }

    int valid = this->InBoundingBox(coords,point);
    valid = valid && this->InCell(coords,point);

    totalVolume += ( valid * this->Volumes.Get(index) );
    return valid;

    }

private:

  int InBoundingBox(CellCoords& coords, dax::Vector3 point) const
  {
    //compute bounds
    dax::Scalar bounds[6] = { coords[0][0],
                              coords[0][0],
                              coords[0][1],
                              coords[0][1],
                              coords[0][2],
                              coords[0][2] };
    for( int i=1; i < 4; ++i )
      {
      bounds[0] = dax::math::Min(coords[i][0], bounds[0]);
      bounds[1] = dax::math::Max(coords[i][0], bounds[1]);

      bounds[2] = dax::math::Min(coords[i][1], bounds[2]);
      bounds[3] = dax::math::Max(coords[i][1], bounds[3]);

      bounds[4] = dax::math::Min(coords[i][2], bounds[4]);
      bounds[5] = dax::math::Max(coords[i][2], bounds[5]);
      }

    int valid = 1;
    // Check if the point is inside the bounding-box
    for( int i=0; i < 3; ++i )
      { valid &= (point[i] >= bounds[i*2]) && (point[i] <= bounds[i*2+1]); }

    return valid;
  }

  //Determines if a point is inside a the cell
  int InCell(CellCoords& coords, dax::Vector3 point) const
  {
    dax::Scalar first = dax::utilities::Determinant(
                  coords[0], coords[1],
                  coords[2], coords[3]);
    dax::Scalar second = dax::utilities::Determinant(
                  point, coords[1],
                  coords[2], coords[3]);

    ( dax::math::SignBit(first) == dax::math::SignBit(second) ) &&
      ( second = dax::utilities::Determinant(
                  coords[0], point,
                  coords[2], coords[3]) );

    ( dax::math::SignBit(first) == dax::math::SignBit(second) ) &&
      ( second = dax::utilities::Determinant(
                  coords[0], coords[1],
                  point, coords[3]) );

    ( dax::math::SignBit(first) == dax::math::SignBit(second) ) &&
      ( second = dax::utilities::Determinant(
                  coords[0], coords[1],
                  coords[2], point) );

    //we care if all the results have the same sign bit
    int valid = dax::math::SignBit(first) == dax::math::SignBit(second);
    return dax::math::Min(valid,1); //make sure we don't have a large valid flag
  }
};

}
}

#endif
