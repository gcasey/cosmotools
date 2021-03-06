#include "LagrangianTesselator.h"

#include "ExtentUtilities.h"
#include "TetrahedronUtilities.h"

// C/C++ includes
#include <cassert> // For assert
#include <cstdlib> // For definition of NULL
#include <cstring> // For memcpy
#include <set>     // For STL set
#include <sstream> // For string streams

#ifdef USEDAX
#include <dax/cont/DeviceAdapter.h>

#include <dax/CellTag.h>
#include <dax/TypeTraits.h>
#include <dax/Extent.h>

#include <dax/cont/ArrayHandle.h>
#include <dax/cont/ArrayHandleConstantValue.h>
#include <dax/cont/GenerateTopology.h>
#include <dax/cont/Scheduler.h>
#include <dax/cont/UniformGrid.h>
#include <dax/cont/UnstructuredGrid.h>

#include <dax/worklet/Tetrahedralize.h>
#endif

namespace cosmologytools
{

struct FaceKey {

  INTEGER ReferenceTetIdx;
  INTEGER ReferenceFaceIdx;
  std::set<INTEGER> FaceIds;

  /**
   * @brief Constructor of a face key
   * @param face the triangular/ face
   */
  FaceKey( INTEGER face[3], INTEGER tetIdx, INTEGER faceIdx )
    {
    this->ReferenceTetIdx  = tetIdx;
    this->ReferenceFaceIdx = faceIdx;

    for( int i=0; i < 3; ++i)
      {
      this->FaceIds.insert( face[i] );
      }
    assert("pre: Face Ids do not correspond to a triangular face" &&
           (this->FaceIds.size()==3));
    }

  /**
   * @brief Gets a key used to uniquely identify this face instance.
   * @return key
   */
  std::string GetKey()
    {
    std::ostringstream oss;
    oss.clear();
    oss.str("");
    std::set<INTEGER>::iterator iter = this->FaceIds.begin();
    for( ; iter != this->FaceIds.end(); ++iter )
      {
      oss << *iter << ".";
      }
    return(oss.str());
    }

};

//=============================================================================

LagrangianTesselator::LagrangianTesselator()
{
  this->Connectivity = NULL;
  this->NumTets      = 0;
  this->Initialize();
}

//-----------------------------------------------------------------------------
LagrangianTesselator::~LagrangianTesselator()
{
  if( this->Connectivity != NULL )
    {
    delete [] this->Connectivity;
    }
  this->Initialize();
}

//-----------------------------------------------------------------------------
void LagrangianTesselator::SetGridExtent(INTEGER extent[6])
{
  memcpy(this->Extent,extent,sizeof(INTEGER)*6);
  INTEGER numVoxels = ExtentUtilities::ComputeNumberOfCells( extent );
  this->NumTets = 5*numVoxels;
  if( this->Connectivity != NULL )
    {
    delete [] this->Connectivity;
    }
  this->Connectivity = new INTEGER[ 4*this->NumTets ];
}

//-----------------------------------------------------------------------------
void LagrangianTesselator::Initialize()
{
  this->Origin[0]  = this->Origin[1]  = this->Origin[2]  = 0.0;
  this->Spacing[0] = this->Spacing[1] = this->Spacing[2] = 0.5;

  for( int i=0; i < 6; ++i )
    {
    this->Extent[i] = 0;
    }

  this->NumTets      = 0;
  this->Connectivity = NULL;
}

//------------------------------------------------------------------------------
void LagrangianTesselator::Tesselate()
{
#ifdef USEDAX
  this->DAXBuildTesselation();
#else
  this->BuildTesselation();
#endif
  this->BuildFaceAdjacency();
}

//------------------------------------------------------------------------------
void LagrangianTesselator::BuildFaceAdjacency()
{
  assert("pre: number of tets" && (this->NumTets > 0) );

  INTEGER tet[4];
  INTEGER face[3];

  for( INTEGER i=0; i < this->NumTets; ++i )
    {
    this->GetTetConnectivity( i, tet );
    for(int fidx=0; fidx < 4; ++fidx )
      {
      TetrahedronUtilities::GetFace(fidx,tet,face);
      FaceKey faceHash(face,i,fidx);

      if(this->FaceAdjacency.find( faceHash.GetKey() ) !=
          this->FaceAdjacency.end() )
        {
        this->FaceAdjacency[ faceHash.GetKey() ].push_back( i );
        }
      else
        {
        std::vector< INTEGER > tets;
        tets.push_back( i );
        this->FaceAdjacency[ faceHash.GetKey() ] = tets;

        std::vector<INTEGER> orientation;
        orientation.push_back( i );
        orientation.push_back( fidx );
        this->FaceToTetOrientation[ faceHash.GetKey() ] = orientation;
        }
      } // END for all tetrahedral faces
    } // END for all tets
}

//------------------------------------------------------------------------------
void LagrangianTesselator::GetBounds(REAL bounds[6], INTEGER fringe)
{
  REAL min[3];
  INTEGER minijk[3];

  REAL max[3];
  INTEGER maxijk[3];

  minijk[0] = IMIN(this->Extent)+fringe;
  minijk[1] = JMIN(this->Extent)+fringe;
  minijk[2] = KMIN(this->Extent)+fringe;
  this->GetPoint(minijk,min);

  maxijk[0] = IMAX(this->Extent)-fringe;
  maxijk[1] = JMAX(this->Extent)-fringe;
  maxijk[2] = KMAX(this->Extent)-fringe;
  this->GetPoint(maxijk,max);

  bounds[0] = min[0];
  bounds[1] = max[0];
  bounds[2] = min[1];
  bounds[3] = max[1];
  bounds[4] = min[2];
  bounds[5] = max[2];
}

//------------------------------------------------------------------------------
void LagrangianTesselator::GetPoint(INTEGER ijk[3], REAL pnt[3])
{
  for( int i=0; i < 3; ++i )
    {
    pnt[ i ] = this->Origin[ i ] + this->Spacing[ i ]*ijk[i];
    } // END for all dimensions
}

#define EXTRACT_LANGRANGIAN_NODE(i,V,T) {                             \
 INTEGER nodeIdx = this->Connectivity[tetIdx*4+i];                    \
 tet[i] = nodeIdx;                                                     \
 INTEGER ijk[3];                                                       \
 ExtentUtilities::GetStructuredCoordinates(nodeIdx,this->Extent,ijk); \
 this->GetPoint( ijk, V);                                              \
}

//------------------------------------------------------------------------------
void LagrangianTesselator::GetLagrangianTet(
    const INTEGER tetIdx,
    REAL v0[3], REAL v1[3], REAL v2[3],REAL v3[3],
    INTEGER tet[4])
{
  assert("pre: tetIdx is out-of-bounds" && (tetIdx < this->NumTets) );

  int i=0;
  EXTRACT_LANGRANGIAN_NODE(i,v0,tet);

  ++i;
  EXTRACT_LANGRANGIAN_NODE(i,v1,tet);

  ++i;
  EXTRACT_LANGRANGIAN_NODE(i,v2,tet);

  ++i;
  EXTRACT_LANGRANGIAN_NODE(i,v3,tet);
}

//------------------------------------------------------------------------------
void LagrangianTesselator::GetTetConnectivity(
        const INTEGER tetIdx, INTEGER tet[4])
{
  assert("pre: tetIdx is out-of-bounds" && (tetIdx < this->NumTets) );

  for( int i=0; i < 4; ++i )
    {
    tet[i] = this->Connectivity[tetIdx*4+i];
    }
}

//------------------------------------------------------------------------------
void LagrangianTesselator::GetFaces(std::vector<INTEGER> &faces)
{
  faces.resize(3*this->FaceAdjacency.size());
  std::map< std::string,std::vector<INTEGER> >::iterator iter =
      this->FaceAdjacency.begin();

  INTEGER tet     = -1;
  INTEGER tetFace = -1;
  INTEGER mytet[4];
  INTEGER myface[3];
  for(int idx=0 ; iter != this->FaceAdjacency.end(); ++iter, ++idx )
    {
    assert(this->FaceToTetOrientation.find(iter->first) !=
           this->FaceToTetOrientation.end());

    std::vector<INTEGER> orientation = this->FaceToTetOrientation[iter->first];
    assert(orientation.size()==2);

    tet     = orientation[ 0 ];
    tetFace = orientation[ 1 ];
    this->GetTetConnectivity(tet,mytet);

    TetrahedronUtilities::GetFace(tetFace,mytet,myface);

    faces[idx*3]   = myface[0];
    faces[idx*3+1] = myface[1];
    faces[idx*3+2] = myface[2];
    } // END for all faces
}

//------------------------------------------------------------------------------
void LagrangianTesselator::GetAdjacentTets(
        INTEGER face[3], std::vector<INTEGER> &tets)
{
  FaceKey faceHash(face,0,0);
  assert("pre: face cannot be found in the FaceAdjacency table!" &&
 (this->FaceAdjacency.find(faceHash.GetKey()) != this->FaceAdjacency.end()));

  tets = this->FaceAdjacency[ faceHash.GetKey() ];

  // sanity check!
  assert("post: face can be shared by at most two tets!" &&
         ((tets.size()==1) || (tets.size()==2)) );
}

//-----------------------------------------------------------------------------
#ifdef USEDAX
void LagrangianTesselator::DAXBuildTesselation()
{
	// STEP 0: Setup  dax Uniform Grid
  const dax::Id3 minE(IMIN(this->Extent), JMIN(this->Extent), KMIN(this->Extent));
  const dax::Id3 maxE(IMAX(this->Extent), JMAX(this->Extent), KMAX(this->Extent));

  dax::cont::UniformGrid<> inputGrid;
  inputGrid.SetExtent(minE,maxE);
  inputGrid.SetOrigin( dax::make_Vector3(this->Origin[0],
                                         this->Origin[1], this->Origin[2]));
  inputGrid.SetSpacing( dax::make_Vector3(this->Spacing[0],
                                          this->Spacing[1], this->Spacing[2]));

  typedef dax::cont::GenerateTopology<dax::worklet::Tetrahedralize,
                                      dax::cont::ArrayHandleConstantValue<dax::Id>
                                      > GenerateT;
  typedef GenerateT::ClassifyResultType  ClassifyResultType;

  dax::cont::Scheduler<> scheduler;

  //for every input cell generate 5 output cells, using a constant
  //value array handle
  ClassifyResultType classification(5,inputGrid.GetNumberOfCells());

  dax::worklet::Tetrahedralize tetWorklet;
  GenerateT generateTets(classification,tetWorklet);

  //don't remove duplicate points.
  generateTets.SetRemoveDuplicatePoints(false);

  //we want to use our Connectivity array as the dax unstructured grids
  //connectivity array. This way dax will directly write into our allocated
  //memory and we won't have to copy anything
  dax::cont::ArrayHandle< dax::Id > connHandle;

  dax::cont::UnstructuredGrid< dax::CellTagTetrahedron > outGrid;
  outGrid.SetCellConnections(connHandle);

  scheduler.Invoke(generateTets,inputGrid,outGrid);

  //this will copy the memory back to the host if we are using cuda, and
  //for OpenMP and TBB it is a no op
  connHandle.CopyInto(this->Connectivity);


}
#endif

//------------------------------------------------------------------------------
#define ADDTET(tetIdx, V0, V1, V2, V3 )   \
    this->Connectivity[tetIdx*4]   = V0;  \
    this->Connectivity[tetIdx*4+1 ]= V1;  \
    this->Connectivity[tetIdx*4+2 ]= V2;  \
    this->Connectivity[tetIdx*4+3 ]= V3;  \
    ++tetIdx;


//-----------------------------------------------------------------------------
void LagrangianTesselator::BuildTesselation()
{
  assert("pre: grid extent is not a 3-D extent!" &&
          ExtentUtilities::Is3DExtent(this->Extent));
  assert("pre: tetrahedral connectivity array is NULL" &&
         (this->Connectivity != NULL) );

  INTEGER V[8];
  INTEGER tetIdx  = 0;
  INTEGER cellIdx = 0;
  INTEGER ijk[3];

  for(I(ijk)=IMIN(this->Extent); I(ijk) < IMAX(this->Extent); ++I(ijk) )
    {
    for(J(ijk)=JMIN(this->Extent); J(ijk) < JMAX(this->Extent); ++J(ijk) )
      {
      for(K(ijk)=KMIN(this->Extent); K(ijk) < KMAX(this->Extent); ++K(ijk) )
        {
        cellIdx = ExtentUtilities::GetLinearIndex(ijk,this->Extent);

        // Nodes on voxel base
        V[0] = ExtentUtilities::GetLinearIndex(
                 I(ijk),J(ijk),K(ijk),this->Extent);
        V[1] = ExtentUtilities::GetLinearIndex(
                 I(ijk)+1,J(ijk),K(ijk),this->Extent);
        V[2] = ExtentUtilities::GetLinearIndex(
                 I(ijk)+1,J(ijk)+1,K(ijk),this->Extent);
        V[3] = ExtentUtilities::GetLinearIndex(
                 I(ijk),J(ijk)+1,K(ijk),this->Extent);

        // Nodes on voxel top
        V[4] = ExtentUtilities::GetLinearIndex(
                 I(ijk),J(ijk),K(ijk)+1,this->Extent);
        V[5] = ExtentUtilities::GetLinearIndex(
                 I(ijk)+1,J(ijk),K(ijk)+1,this->Extent);
        V[6] = ExtentUtilities::GetLinearIndex(
                 I(ijk)+1,J(ijk)+1,K(ijk)+1,this->Extent);
        V[7] = ExtentUtilities::GetLinearIndex(
                 I(ijk),J(ijk)+1,K(ijk)+1,this->Extent);

        // Set tetrahedral connectivity
        if( cellIdx % 2 )
          {
          ADDTET(tetIdx,V[1],V[3],V[6],V[2]);
          ADDTET(tetIdx,V[0],V[1],V[3],V[4]);
          ADDTET(tetIdx,V[1],V[4],V[5],V[6]);
          ADDTET(tetIdx,V[3],V[6],V[7],V[4]);
          ADDTET(tetIdx,V[1],V[4],V[6],V[3]);
          }
        else
          {
          ADDTET(tetIdx,V[2],V[1],V[5],V[0]);
          ADDTET(tetIdx,V[0],V[2],V[3],V[7]);
          ADDTET(tetIdx,V[2],V[5],V[6],V[7]);
          ADDTET(tetIdx,V[0],V[7],V[4],V[5]);
          ADDTET(tetIdx,V[0],V[2],V[7],V[5]);
          }

        } // END for all k
      } // END for all j
    } // END for all i
}

} /* namespace cosmologytools */
