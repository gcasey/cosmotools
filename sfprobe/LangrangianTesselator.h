/**
 * @class LagrangianTesselator
 * @brief A class that implements template-based tesselation of a user-supplied
 * extent.
 */
#ifndef LangrangianTesselator_H_
#define LangrangianTesselator_H_

#include "CosmologyToolsMacros.h"

namespace cosmologytools
{

class LangrangianTesselator
{
public:
  LangrangianTesselator();
  virtual ~LangrangianTesselator();

  /**
   * @brief Set/Get the grid origin
   */
  SetVector3Macro(Origin,REAL);
  GetMacro(Origin, const REAL*);

  /**
   * @brief Set/Get the grid spacing
   */
  SetVector3Macro(Spacing,REAL);
  GetMacro(Spacing,const REAL*);

  /**
   * @brief Sets the virtual grid extent for this instance.
   * @param extent the node extent for this StructureFormationProbe instance.
   * @note extent is a flat array that stores[imin imax jmin jmax kmin kmax]
   * @pre extent is expected to be a 3-D extent
   */
  void SetGridExtent(int extent[6]);
  GetMacro(Extent,const INTEGER*);

  /**
   * @brief Returns the total number of tetrahedra generate
   */
  GetMacro(NumTets,INTEGER);

  /**
   * @brief Returns the vertex ids w.r.t. the given extent that compose the tet
   * of the given tetrahedron index.
   * @param tetIdx the index of the tetrahedron
   * @param tet the list of vertex IDs that make up the tet (out)
   */
  void GetTetConnectivity(const INTEGER tetIdx, INTEGER tet[4]);

  /**
   * @brief Returns the tetrahedron in langrangian space
   * @param tetIdx the tetrahedral index
   * @param v0 the langrangian coordinates of the 1st vertex
   * @param v1 the langrangian coordinates of the 2nd vertex
   * @param v2 the langrangian coordinates of the 3rd vertex
   * @param v3 the langrangian coordinates of the 4th vertex
   * @param tet the tetrahedral connectivity
   */
  void GetLangrangianTet(
      const INTEGER tetIdx,
      REAL v0[3], REAL v1[3], REAL v2[3],REAL v3[3],
      INTEGER tet[4]
      );


protected:
  REAL Origin[3];         // The origin of the grid
  REAL Spacing[3];        // The spacing of the grid

  INTEGER Extent[6];      // The grid node extent (supplied by the user)

  INTEGER NumTets;        // Total number of tets
  INTEGER *Connectivity;  // tet connectivity list


  /**
   * @brief Get the coordinates of the point in langrangian space
   * @param ijk the structured coordinates of the point
   * @param pnt the cartesian coordinates of the point
   */
  void GetPoint(INTEGER ijk[3], REAL pnt[3]);

  /**
   * @brief Constructs tesselation in langrangian space.
   */
  void BuildTesselation();

  /**
   * @brief Initialize
   */
  void Initialize();

private:
  DISABLE_COPY_AND_ASSIGNMENT(LangrangianTesselator);


};

} /* namespace cosmologytools */
#endif /* LangrangianTesselator_H_ */
