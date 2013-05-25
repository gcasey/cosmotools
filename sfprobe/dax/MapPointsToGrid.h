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
#ifndef __MapPointsToGrid_worklet_
#define __MapPointsToGrid_worklet_

#include <dax/Extent.h>
#include <dax/math/Precision.h>
#include <dax/exec/WorkletMapField.h>

namespace dax {
namespace worklet {

class MapPointsToGrid : public dax::exec::WorkletMapField
{
public:
  //origin, spacing, dims
  typedef void ControlSignature(Field(In), Field(In), Field(In), Field(In), Field(Out));
  typedef _5 ExecutionSignature(_1, _2, _3, _4, WorkId);

  DAX_EXEC_EXPORT
  dax::Tuple<dax::Id,2> operator()(const dax::Vector3& coordinate,
                     const dax::Vector3& origin,
                     const dax::Vector3& spacing,
                     const dax::Extent3& extent,
                     dax::Id id) const
  {
  dax::Tuple<dax::Id, 2> result;
  result[0] = -1;
  result[1] = id;

  //for each point determine which bucket in the virtual grid it is inside.
  const dax::Id3 dist( coordinate[0] - origin[0],
                       coordinate[1] - origin[1],
                       coordinate[2] - origin[2]);


  const dax::Id3 ijk( dax::math::Floor(dist[0] / spacing[0]),
                      dax::math::Floor(dist[1] / spacing[1]),
                      dax::math::Floor(dist[2] / spacing[2]));

  if( (ijk[0] < 0 || ijk[1] < 0 || ijk[2] < 0) ||
    (ijk[0] > extent.Max[0] || ijk[1] > extent.Max[1] || ijk[2] > extent.Max[2] ) )
    {
    return result;
    }

  //confirm we have a ijk that is all positive
  dax::Id li = ijk[0]-extent.Min[0];
  dax::Id lj = ijk[1]-extent.Min[1];
  dax::Id lk = ijk[2]-extent.Min[2];

  // Get Extent dimensions
  dax::Id N2  = extent.Max[1]-extent.Min[1]+1;
  dax::Id N1  = extent.Max[0]-extent.Min[0]+1;

  // Return linear index
  result[0] = ( (lk*N2+lj)*N1+li );
  return result;
  }
};

}
}


#endif
