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
#ifndef __MapTetsToPoints_h
#define __MapTetsToPoints_h

#include "TetStructure.h"
#include "ConstArray.h"

#include <dax/Extent.h>
#include <dax/math/Compare.h>
#include <dax/exec/WorkletMapField.h>

namespace dax {
namespace worklet {

class MapTetsToPoints : public dax::exec::WorkletMapField
{
public:
  //origin, spacing, dims
  typedef void ControlSignature(Field(In), ExecObject(In), ExecObject(In), Field(Out), Field(Out));
  typedef void ExecutionSignature(_1, _2, _3, _4, _5);

  DAX_EXEC_EXPORT
  void operator()(const dax::Vector3& coordinate,
                  const dax::exec::TetStructure& structure,
                  const dax::exec::ConstArray<dax::Id> tetIdsToSearch,
                  dax::Id& stream,
                  dax::Scalar& rho) const
  {
  //for each point determine which cells in the tet structure are inside it
  const dax::Id numTets = tetIdsToSearch.GetNumberOfValues();

  dax::Id numValidCells = 0;
  dax::Scalar volumeSum = 0;

  for( dax::Id i = 0; i != numTets; ++i)
    {
    dax::Id tetId  = tetIdsToSearch.Get(i);
    numValidCells += structure.ProbePoint(tetId, coordinate, volumeSum);
    }

  //write back to memory
  stream = numValidCells;
  if(numValidCells > 0)
    {
    rho = ( volumeSum / static_cast<dax::Scalar>(numValidCells) );
    }
  }
};

}
}

#endif
