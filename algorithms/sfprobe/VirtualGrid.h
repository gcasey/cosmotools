/**
 * @class VirtualGrid
 * @brief Simple implementation of virtual grid used to minimize the search
 * space for probing a mesh.
 */
#ifndef VIRTUALGRID_H_
#define VIRTUALGRID_H_

#include "CosmoToolsMacros.h"

// CosmologyTools includes
#include "ExtentUtilities.h"

// C/C++ includes
#include <vector>
#include <set>
#include <string>

namespace cosmologytools {

class SimpleMesh;

class VirtualGrid
{
public:
  VirtualGrid();
  virtual ~VirtualGrid();

  SetVector3Macro(Dimensions,INTEGER);
  GetMacro(Dimensions,const INTEGER*)
  GetMacro(Spacing,const REAL*);
  GetMacro(Bounds, const REAL*);

  int GetNumberOfBuckets() const { return Buckets.size(); }

  std::set<INTEGER> const& GetBucket( int i) { return Buckets[i]; }

  /**
   * @brief Registers the mesh
   * @param M the user-supplied mesh
   */
  void RegisterMesh(SimpleMesh &M);

  /**
   * @brief Return the candiate cells for the given query point
   * @param pnt the query point
   * @param cells the list of candidate cells that contain the point.
   */
  void GetCandidateCellsForPoint(
        REAL pnt[3],std::vector<INTEGER> &cells);

  /**
   * @brief Returns this virtual grid as a legacy VTK grid.
   * @return s a string representation of this instance of virtual grid.
   */
  std::string ToLegacyVtkString();


protected:
  INTEGER Dimensions[3];

  REAL Bounds[6];
  REAL Spacing[3];
  std::vector< std::set<INTEGER> > Buckets;

  void GetCellExtent(INTEGER ext[6])
    {
    IMIN(ext) = JMIN(ext) = KMIN(ext) = 0;
    IMAX(ext) = this->Dimensions[0]-1;
    JMAX(ext) = this->Dimensions[1]-1;
    KMAX(ext) = this->Dimensions[2]-1;
    }

  /**
   * @brief Gets the (local) extent corresponding to this virtual grid instance.
   * @param ext the extent of the virtual grid (out)
   */
  void GetExtent(INTEGER ext[6])
    {
    IMIN(ext) = JMIN(ext) = KMIN(ext) = 0;
    IMAX(ext) = this->Dimensions[0];
    JMAX(ext) = this->Dimensions[1];
    KMAX(ext) = this->Dimensions[2];
    };

  /**
   * @brief Injects a cell to a bucket in the virtual grid.
   * @param cellIdx the index of the cell to inject
   * @param M reference to the mesh that has the cell
   * @return status true iff the cell is succesfully injected in to a bucket of
   * the virtual grid, else, false.
   */
  bool InjectCell(INTEGER cellIdx, SimpleMesh &M);

  /**
   * @brief Given the cartesian coordinates of the point, this method finds the
   * bucket in the virtual grid that contains the point.
   * @param pnt the cartesian coordinates of the query point (in)
   * @param ijk the structured coordinates of the bucket (out)
   * @return status true if the point is inside the virtual grid, else false.
   */
  bool FindBucket(REAL pnt[3], INTEGER ijk[3]);

  /**
   * @brief Checks if the given IJK coordinates are within the cell extent of
   * this virtual grid instance.
   * @param ijk the ijk coordinates (in)
   * @return status true if ijk is within the cell extent, else false.
   */
  bool WithinCellExtent(INTEGER ijk[3])
    {
    INTEGER ext[6];
    this->GetExtent( ext );
    for( int i=0; i < 3; ++i)
      {
      // NOTE: ext[i*2+1]-1 is used here b/c we are checking if the ijk is
      // within the cell-extent of the virtual grid.
      if( (ijk[i] < ext[i*2]) || (ijk[i] > ext[i*2+1]-1) )
        {
        return false;
        }
      } // END for all dimensions
    return true;
    };

private:
  DISABLE_COPY_AND_ASSIGNMENT(VirtualGrid);
};

} /* namespace cosmologytools */
#endif /* VIRTUALGRID_H_ */
