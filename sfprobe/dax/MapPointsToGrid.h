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

#include <dax/exec/WorkletMapField.h>
#include <dax/Extent.h>
#include "dax/Types.h"

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
                     const dax::Extent3& extents,
                     dax::Id id) const
  {
  //for each point determine which cell in the virtual grid it is inside.
  const dax::Vector3 vijk( (coordinate - origin) / spacing );
  const dax::Id3 ijk( vijk[0], vijk[1], vijk[2] );

  dax::Tuple<dax::Id, 2> result;
  result[0] = dax::index3ToFlatIndexCell(ijk,extents);
  result[1] = id;
  return result;
  }
};

}
}

#endif
