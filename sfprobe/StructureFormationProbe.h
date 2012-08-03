/**
 * @class StructureFormationProbe
 * @brief A class the implements a tesselation-based approach for finding
 * streams and caustic surfaces in N-Body cosmology simulations.
 */
#ifndef STRUCTUREFORMATIONPROBE_H_
#define STRUCTUREFORMATIONPROBE_H_

#include "CosmologyToolsMacros.h"
#include "LangrangianTesselator.h"
#include "SimpleMesh.h"

// C/C++ includes
#include <vector>
#include <map>

namespace cosmologytools {


class StructureFormationProbe
{
public:
  StructureFormationProbe();
  virtual ~StructureFormationProbe();

  /**
   * @brief Get/Set macro for fringe
   */
  GetNSetMacro(Fringe,INTEGER);

  /**
   * @brief Sets the particle (positions) and corresponding global Ids.
   * @param particles the particles position vector strided by 3.
   * @param GlobalIds list of the global Ids corresponding to each particle.
   * @param N the total number of particles.
   * @note the particle global IDs must start numbering from 0.
   */
  void SetParticles(REAL *particles, INTEGER *GlobalIds, INTEGER N);

  /**
   * @brief Constructs the langrangian mesh with the given origin, spacing
   * and extent.
   */
  void BuildLangrangianMesh(
          REAL origin[3], REAL spacing[3], INTEGER extent[6]);


  /**
   * @brief Returns the langrangian tesselator instance that is associated with
   * this instance of StructureFormationProbe.
   * @return langrange
   */
  LangrangianTesselator* GetLangrangeTesselator();

  /**
   * @brief Builds the euler mesh by mapping the tet connectivity of the
   * langrangian mesh onto the given particles.
   * @pre A langragian tesselation must have been first constructed.
   */
  void BuildEulerMesh();

  /**
   * @brief Returns the nodes, tet connectivity of the computed euler mesh as
   * well as the volume for each element in the mesh.
   * @param nodes the euler mesh nodes (out)
   * @param tets the euler mesh tetrahedral connectivity (out)
   * @param volumes the volumes of each cell in the euler mesh (out)
   */
  void GetEulerMesh(
        std::vector<REAL> &nodes,
        std::vector<INTEGER> &tets,
        std::vector<REAL> &volumes);

  /**
   * @brief Loops through all mesh faces and checks if the face corresponds to
   * a face on a caustic surface.
   * @param nodes
   * @param triangles
   */
  void ExtractCausticSurfaces(
      std::vector<REAL> &nodes, std::vector<INTEGER> &triangles);

  /**
   * @brief Find the number of streams at the given location.
   * @param pnt the xyz coordinates of the probe point.
   * @return N the number of streams at the user-supplied location.
   * @note The number of streams is equal to the number of tetrahedra that
   * contain the given point.
   */
  INTEGER GetNumberOfStreams(REAL pnt[3]);

protected:
  REAL *Particles;      // pointer to user-supplied array of particles
  INTEGER *GlobalIds;   // pointer to user-supplied array of global IDs
  INTEGER NumParticles; // the total number of particles (user-supplied)
  INTEGER Fringe;       // user-supplied parameter to reject tets that are
                        // fringe distance from the boundary. The fringe is
                        // used to deal with artifacts on the periodic boundary.

  /**
   * @brief Given the IDs of two tetrahedra, this methods returns true if their
   * respective volumes have opposite signs.
   * @param tetIdx1 the ID of the first tet
   * @param tetIdx2 the ID of the second tet
   * @return status true if the volumes of the two tets have opposite volume,
   * else false.
   */
  bool VolumesHaveOppositeSigns(INTEGER tetIdx1, INTEGER tetIdx2);

  /**
   * @brief Given the tet connectivity in langrangian space, this method returns
   * the xyz coordinates of the tet nodes in euler-space corresponding to the
   * user-supplied particles vector. If the tet is mapped successfully this
   * method returns true. Otherwise, false is returned.
   * @param tet the tetrahedral connectivity
   * @param nodes flat array where the nodes will be stored.
   * @return status true if successful, else false.
   */
  bool MapTetToEulerSpace(INTEGER tet[4], REAL nodes[12]);

  /**
   * @brief Checks if the node is within an interior box constructed based on the
   * langrangian grid box extent domain.
   * @param nodeIdx
   * @return
   */
  bool IsNodeWithinFringeBounds(INTEGER nodeIdx);

  /**
   * @brief This method checks if the given node is within a fringe from a
   * periodic boundary.
   * @param nodeIdx
   * @return
   */
  bool IsNodeWithinPeriodicBoundaryFringe( INTEGER nodeIdx );

  // Mapping of global particle IDs to the particle storage locations in the
  // user-supplied array.
  std::map< INTEGER, INTEGER > Global2PositionMap;

  // Mapping of Langrange mesh entities (nodes/cells) to the euler mesh
  std::map<INTEGER,INTEGER> LangrangeNode2EulerNode;
  std::map< INTEGER, INTEGER > LangrangeTet2EulerTet;

  LangrangianTesselator *Langrange; // langrangian tesselator (computed)

  SimpleMesh EulerMesh;
  std::vector< REAL > Volumes;

private:
  DISABLE_COPY_AND_ASSIGNMENT(StructureFormationProbe);
};

} /* namespace hacctools */
#endif /* STRUCTUREFORMATIONPROBE_H_ */
