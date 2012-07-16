/**
 * @class ExtentUtilities
 * @brief A singleton class that implements fuctionality for 3-D extent
 * calculations.
 */
#ifndef EXTENTUTILITIES_H_
#define EXTENTUTILITIES_H_

// Some usefull extent macros
#define IMIN(ext) ext[0]
#define IMAX(ext) ext[1]
#define JMIN(ext) ext[2]
#define JMAX(ext) ext[3]
#define KMIN(ext) ext[4]
#define KMAX(ext) ext[5]

// Some useful IJK macros
#define I(ijk) ijk[0]
#define J(ijk) ijk[1]
#define K(ijk) ijk[2]

#include "CosmologyToolsMacros.h"

namespace cosmologytools
{

class ExtentUtilities
{
public:
  ExtentUtilities();
  virtual ~ExtentUtilities();

  /**
   * @brief Checks if this extent is a 3-D extent.
   * @param ext the extent to check
   * @return true if the extent is a 3-D extent else false.
   */
  static bool Is3DExtent( INTEGER ext[6] );

  /**
   * @brief Given structured coordinates of a point, this method computes the
   * linear index w.r.t. the given extent.
   * @param ijk the structured coordinates of the point
   * @param ext the extent that the linear index is calculated w.r.t. to
   * @return idx the linear index of the node w.r.t. the supplied extent.
   */
  static INTEGER GetLinearIndex(INTEGER ijk[3], INTEGER ext[6]);
  static INTEGER GetLinearIndex(
      const INTEGER i, const INTEGER j, const INTEGER k, INTEGER ext[6]);

  /**
   * @brief Computes the structured coordinates of the point with the given
   * linear index w.r.t. to the given grid extent.
   * @param idx the linear index of the point
   * @param ext the extent of the grid
   * @param ijk the corresponding IJK coordinates
   */
  static void GetStructuredCoordinates(
      INTEGER idx, INTEGER ext[6], INTEGER ijk[3]);

  /**
   * @brief Computes the number of nodes in the given extent.
   * @param ext the user-supplied extent.
   * @return N the number of nodes in the grid corresponding to the given
   * extent.
   * @pre this->Is3DExtent( ext ) == true
   * @post N > 0
   */
  static INTEGER ComputeNumberOfNodes( INTEGER ext[6] );

  /**
   * @brief Computes the total number of cells in the given extent
   * @param ext the user-supplied extent
   * @return N the number of cells in the grid corresponding to the
   * given extent.
   * @pre this->Is3DExtent( ext ) == true
   * @post N > 0
   */
  static INTEGER ComputeNumberOfCells( INTEGER ext[6] );

private:
  DISABLE_COPY_AND_ASSIGNMENT(ExtentUtilities);
};

} /* namespace cosmologytools */
#endif /* EXTENTUTILITIES_H_ */
