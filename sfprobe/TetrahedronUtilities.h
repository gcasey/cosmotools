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
