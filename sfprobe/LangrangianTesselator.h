/**
 * @class LagrangianTesselator
 * @brief A class that implements template-based tesselation of a user-supplied
 * extent.
 */
#ifndef LangrangianTesselator_H_
#define LangrangianTesselator_H_

#include "CosmologyToolsMacros.h"

// C/C++ includes
#include <map>    // For STL map
#include <vector> // For STL vector
#include <string> // For C++ string

namespace cosmologytools
{

// Forward Declarations
struct FaceKey;

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
   * @brief Returns the bounds of langrangian domain
   * @param bounds user-supplied array where the bounds will be stored
   * @param fringe optional parameter to get interior fringe bounds
   * @note the bounds are stored [xmin xmax ymin ymax zmin zmax]
   */
void GetBounds(REAL bounds[6], INTEGER fringe=0);

  /**
   * @brief Constructs the tesselation.
   */
  void Tesselate();

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

  /**
   * @brief Return the triangular faces.
   * @param faces the output faces vector (strided by 3)
   */
  void GetFaces(std::vector<INTEGER> &faces);

  /**
   * @brief Returns the adjaces tets of the given face.
   * @param face the Ids that make up the face
   * @param tets the tet Ids that share the commong face
   * @post tets.size()==1 || tets.size()==2
   * @note if tets.size()==1, the tet is on a boundary.
   */
  void GetAdjacentTets(INTEGER face[3], std::vector<INTEGER> &tets);

protected:
  REAL Origin[3];         // The origin of the grid
  REAL Spacing[3];        // The spacing of the grid

  INTEGER Extent[6];      // The grid node extent (supplied by the user)

  INTEGER NumTets;        // Total number of tets
  INTEGER *Connectivity;  // tet connectivity list

  // Data-structure to store the face adjacency information
  std::map< std::string,std::vector<INTEGER> > FaceAdjacency;
  std::map< std::string,std::vector<INTEGER> > FaceToTetOrientation;

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
   * @brief Constructs the face adjacency information.
   */
  void BuildFaceAdjacency();

  /**
   * @brief Initialize
   */
  void Initialize();

private:
  DISABLE_COPY_AND_ASSIGNMENT(LangrangianTesselator);
};

} /* namespace cosmologytools */
#endif /* LangrangianTesselator_H_ */
