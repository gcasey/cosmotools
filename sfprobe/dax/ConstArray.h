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
#ifndef __ConstArray_h
#define __ConstArray_h

#include <dax/cont/ArrayHandle.h>
#include <dax/exec/ExecutionObjectBase.h>
#include <vector>

namespace dax {
namespace exec {

template <typename ValueType>
class ConstArray : public dax::exec::ExecutionObjectBase
{
  typedef typename dax::cont::ArrayHandle< ValueType >::PortalConstExecution PortalType;

  PortalType Portal;
public:
  DAX_CONT_EXPORT ConstArray( dax::cont::ArrayHandle< ValueType > arrayHandle ):
  Portal(arrayHandle.PrepareForInput())
  {
    //by constructing a ConstVector you transfer all the array information
    //from the host to the device machine
  }

  DAX_EXEC_EXPORT dax::Id GetNumberOfValues( ) const
    {
    return this->Portal.GetNumberOfValues();
    }

  DAX_EXEC_EXPORT ValueType Get( dax::Id index ) const
    {
    return this->Portal.Get(index);
    }
};

}
}

#endif
