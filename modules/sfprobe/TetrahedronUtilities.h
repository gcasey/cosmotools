/**
 * @class TetrahedronUtilities
 * @brief A singleton class that provides operators for tetrahedra.
 */
#ifndef TETRAHEDRONUTILITIES_H_
#define TETRAHEDRONUTILITIES_H_

#include "CosmologyToolsMacros.h"

namespace cosmologytools
{

class TetrahedronUtilities
{
public:
  TetrahedronUtilities();
  virtual ~TetrahedronUtilities();

  /**
   * @brief Checks if the point is within the bounding of the tet.
   * @param pnt cartesian coordinates of the query point
   * @param V0 cartesian coordinates of base vector
   * @param V1 cartesian coordinates of base vector
   * @param V2 cartesian coordinates of base vector
   * @param V3 cartesian coordinates at the apex
   * @return status true if the point is within the bounding of the tet
   * @note This method is meant as a fast check to be used before calling
   * HasPoint.
   */
  static bool PointInTetBoundingBox(
      REAL pnt[3],REAL V0[3],REAL V1[3], REAL V2[3], REAL V3[3]);

  /**
   * @brief Checks if the point is contained in the tet with the supplied coords
   * @param pnt cartesian coordinates of the query point
   * @param V0 cartesian coordinates of base vector
   * @param V1 cartesian coordinates of base vector
   * @param V2 cartesian coordinates of base vector
   * @param V3 cartesian coordinates at the apex
   * @return status true if the point is inside the tet, else, false.
   */
  static bool HasPoint(
      REAL pnt[3],REAL V0[3],REAL V1[3], REAL V2[3], REAL V3[3]);

  /**
   * @brief Computes the volume of the tet given cartesian coordinates.
   * @param v0 cartesian coordinates of base vector
   * @param v1 cartesian coordinates of base vector
   * @param v2 cartesian coordinates of base vector
   * @param v3 cartesian coordinates at the apex
   * @return V the volume of the given tet
   */
  static REAL ComputeVolume(
      REAL v0[3], REAL v1[3], REAL v2[3], REAL v3[3]);

  /**
   * @brief Returns the IDs of the requested face
   * @param idx the index of the requested face
   * @param tet the connectivity of the tetrahedron
   * @param face user-supplied storage for the face IDs (out)
   * @pre idx >= 0 && idx < 4
   */
  static void GetFace( const int idx,INTEGER tet[4], INTEGER face[3]);

protected:

  /**
   * @brief Checks if the sign of two numbers a,b is the same
   * @param a user-supplied number to check
   * @param b user-supplied number to check
   * @return status true iff a,b have the same sign, else false.
   */
  static bool SameSign(const REAL &a,const REAL &b)
    { return((a*b)>0.0f)?true:false;};

  /**
   * @brief Computes the determinant of the given 2x2 matrix
   * @note the matrix is given in the following form: <br/>
   * |a c| <br/>
   * |b d| <br/>
   * @return d the determinant of the 2x2 matrix
   */
  static REAL Determinant2x2(REAL a, REAL b, REAL c, REAL d)
    {return(a*d-b*c);};

  /**
   * @brief Computes the determinant of the given 3x3 matrix
   * @note the matrix is given in the following form: <br/>
   * |a1 b1 c1| <br/>
   * |a2 b2 c2| <br/>
   * |a3 b3 c3| <br/>
   * @return d the determinant of the 3x3 matrix
   */
  static REAL Determinant3x3(
      REAL a1, REAL a2, REAL a3,
      REAL b1, REAL b2, REAL b3,
      REAL c1, REAL c2, REAL c3)
      {
        REAL det = ( (a1*TetrahedronUtilities::Determinant2x2(b2,b3,c2,c3) )-
                     (b1*TetrahedronUtilities::Determinant2x2(a2,a3,c2,c3) )+
                     (c1*TetrahedronUtilities::Determinant2x2(a2,a3,b2,b3) )
                     );
        return( det );
      }


  /**
   * @brief Computes the determinant of a 4x4 matrix
   * @note the matrix is given in the following form: <br/>
   * |a1 b1 c1 d1| <br/>
   * |a2 b2 c2 d2| <br/>
   * |a3 b3 c3 d3| <br/>
   * |a4 b4 c4 d4| <br/>
   * @return d the determinannt of the martix
   */
  static REAL Determinant(
      REAL a1, REAL a2, REAL a3, REAL a4,
      REAL b1, REAL b2, REAL b3, REAL b4,
      REAL c1, REAL c2, REAL c3, REAL c4,
      REAL d1, REAL d2, REAL d3, REAL d4)
      {
        REAL det = (
         (a1*TetrahedronUtilities::Determinant3x3(b2,b3,b4,c2,c3,c4,d2,d3,d4))-
         (b1*TetrahedronUtilities::Determinant3x3(a2,a3,a4,c2,c3,c4,d2,d3,d4))+
         (c1*TetrahedronUtilities::Determinant3x3(a2,a3,a4,b2,b3,b4,d2,d3,d4))-
         (d1*TetrahedronUtilities::Determinant3x3(a2,a3,a4,b2,b3,b4,c2,c3,c4))
         );
        return( det );
      }

  /**
   * @brief Computes the dot-product of two vector x,y
   * @param x the x-vector
   * @param y the y-vector
   * @return r the dot product
   */
  static REAL DotProduct(REAL x[3],REAL y[3]) {
    return(x[0] * y[0] + x[1] * y[1] + x[2] * y[2] ); };

  /**
   * @brief Computes the cross-product of two vectors, Z = X (x) Y
   * @param x the x-vector
   * @param y the y-vector
   * @param z the resulting cross-product (out)
   */
  static void CrossProduct(REAL x[3],REAL y[3],REAL z[3])
    {
    z[0] = x[1] * y[2] - x[2] * y[1];
    z[1] = x[2] * y[0] - x[0] * y[2];
    z[2] = x[0] * y[1] - x[1] * y[0];
    }

private:
  DISABLE_COPY_AND_ASSIGNMENT(TetrahedronUtilities);
};

} /* namespace cosmologytools */
#endif /* TETRAHEDRONUTILITIES_H_ */
