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
#ifndef __dax_utilities_functions_h
#define __dax_utilities_functions_h

namespace dax {
namespace utilities {
  /**
   * @brief Computes the determinant of the given 2x2 matrix
   * @note the matrix is given in the following form: <br/>
   * |a c| <br/>
   * |b d| <br/>
   * @return d the determinant of the 2x2 matrix
   */
  DAX_EXEC_EXPORT dax::Scalar Determinant2x2(dax::Scalar a, dax::Scalar b,
                                    dax::Scalar c, dax::Scalar d)
    {return(a*d-b*c);};

  /**
   * @brief Computes the determinant of the given 3x3 matrix
   * @note the matrix is given in the following form: <br/>
   * |a1 b1 c1| <br/>
   * |a2 b2 c2| <br/>
   * |a3 b3 c3| <br/>
   * @return d the determinant of the 3x3 matrix
   */
  DAX_EXEC_EXPORT dax::Scalar Determinant3x3(
      dax::Scalar a1, dax::Scalar a2, dax::Scalar a3,
      dax::Scalar b1, dax::Scalar b2, dax::Scalar b3,
      dax::Scalar c1, dax::Scalar c2, dax::Scalar c3)
      {
        return ( (a1*Determinant2x2(b2,b3,c2,c3) )-
                 (b1*Determinant2x2(a2,a3,c2,c3) )+
                 (c1*Determinant2x2(a2,a3,b2,b3) )
               );
      }


  /**
   * @brief Computes the determinant of a 4x4 matrix
      where all d values are set to 1.of
   * @return d the determinannt of the martix
   */
  DAX_EXEC_EXPORT dax::Scalar Determinant(
      const dax::Vector3& a,
      const dax::Vector3& b,
      const dax::Vector3& c,
      const dax::Vector3& d)
      {
        return (
         (a[0]*Determinant3x3(b[1],c[1],d[1],b[2],c[2],d[2],1.0,1.0,1.0))-
         (a[1]*Determinant3x3(b[0],c[0],d[0],b[2],c[2],d[2],1.0,1.0,1.0))+
         (a[2]*Determinant3x3(b[0],c[0],d[0],b[1],c[1],d[1],1.0,1.0,1.0))-
         (1.0*Determinant3x3(b[0],c[0],d[0],b[1],c[1],d[1],b[2],c[2],d[2]))
         );
      }

} }
#endif
