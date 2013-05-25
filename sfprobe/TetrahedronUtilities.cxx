#include "TetrahedronUtilities.h"

// C/C++ includes
#include <iostream> // For cout/cerr
#include <cassert>  // For assert

namespace cosmologytools
{

TetrahedronUtilities::TetrahedronUtilities()
{
  // TODO Auto-generated constructor stub

}

//-----------------------------------------------------------------------------
TetrahedronUtilities::~TetrahedronUtilities()
{
  // TODO Auto-generated destructor stub
}

//-----------------------------------------------------------------------------
#define CHECK_POINT(V,bounds) {     \
    for(int i=0; i < 3; ++i) {      \
      if(V[i] < bounds[i*2])         \
        bounds[i*2]=V[i];            \
      else if(V[i] > bounds[i*2+1]) \
        bounds[i*2+1]=V[i];          \
    } /* END for all dimensions */   \
}

bool TetrahedronUtilities::PointInTetBoundingBox(
    REAL pnt[3],REAL V0[3],REAL V1[3], REAL V2[3], REAL V3[3])
{
    // STEP 0: Initialize bounds
    REAL bounds[6];
    for( int i=0; i < 3; ++i )
      {
      bounds[i*2]=bounds[i*2+1]=V0[i];
      } // END for all dimensions

    // STEP 1: Compute bounding box
    CHECK_POINT(V1,bounds);
    CHECK_POINT(V2,bounds);
    CHECK_POINT(V3,bounds);

    // STEP 2: Check if the point is inside the bounding-box
    for( int i=0; i < 3; ++i )
      {
      if( (pnt[i] < bounds[i*2]) || (pnt[i] > bounds[i*2+1]) )
        return false;
      } // END for all dimensions

    return true;
}

//-----------------------------------------------------------------------------
bool TetrahedronUtilities::HasPoint(
    REAL pnt[3],REAL V0[3],REAL V1[3], REAL V2[3], REAL V3[3])
{

  REAL d0 = TetrahedronUtilities::Determinant(
                V0[0],V1[0],V2[0],V3[0],
                V0[1],V1[1],V2[1],V3[1],
                V0[2],V1[2],V2[2],V3[2],
                1.,1.,1.,1.
                );

  REAL d1 = TetrahedronUtilities::Determinant(
              pnt[0],V1[0],V2[0],V3[0],
              pnt[1],V1[1],V2[1],V3[1],
              pnt[2],V1[2],V2[2],V3[2],
              1.,1.,1.,1.
              );
  if( !TetrahedronUtilities::SameSign(d0,d1) )
    {
    return false;
    }

  REAL d2 = TetrahedronUtilities::Determinant(
              V0[0],pnt[0],V2[0],V3[0],
              V0[1],pnt[1],V2[1],V3[1],
              V0[2],pnt[2],V2[2],V3[2],
              1.,1.,1.,1.
              );
  if( !TetrahedronUtilities::SameSign(d1,d2))
    {
    return false;
    }

  REAL d3 = TetrahedronUtilities::Determinant(
              V0[0],V1[0],pnt[0],V3[0],
              V0[1],V1[1],pnt[1],V3[1],
              V0[2],V1[2],pnt[2],V3[2],
              1.,1.,1.,1.
              );
  if( !TetrahedronUtilities::SameSign(d2,d3) )
    {
    return false;
    }

  REAL d4 = TetrahedronUtilities::Determinant(
              V0[0],V1[0],V2[0],pnt[0],
              V0[1],V1[1],V2[1],pnt[1],
              V0[2],V1[2],V2[2],pnt[2],
              1.,1.,1.,1.
              );
  if( !TetrahedronUtilities::SameSign(d3,d4))
    {
    return false;
    }

  // All 5 determinants have the same sign
  return true;
}

//-----------------------------------------------------------------------------
REAL TetrahedronUtilities::ComputeVolume(
      REAL v0[3], REAL v1[3], REAL v2[3], REAL v3[3] )
{

  REAL vol = 0.;
  REAL A[3], B[3], C[3], BxC[3];

  for( int i=0; i < 3; ++i)
    {
    // vector pointing from 3 to 0 (apex)
    A[ i ] = v3[ i ] - v0[ i ];

    // vector pointing from 0 to 1
    B[ i ] = v1[ i ] - v0[ i ];

    // vector pointing from 0 to 2
    C[ i ] = v2[ i ] - v0[ i ];
    } // END for all dimensions

  TetrahedronUtilities::CrossProduct(B,C,BxC);
  vol = 0.1666666667 * TetrahedronUtilities::DotProduct(A,BxC);
  return( vol );
}

//-----------------------------------------------------------------------------
void TetrahedronUtilities::GetFace(
        const int idx, INTEGER tet[4], INTEGER face[3])
{
  assert("pre: face index is out-of bounds" &&
          (idx >= 0) && (idx < 4) );

  switch( idx )
    {
    case 0:
      face[0] = tet[0];
      face[1] = tet[3];
      face[2] = tet[1];
      break;
    case 1:
      face[0] = tet[1];
      face[1] = tet[3];
      face[2] = tet[2];
      break;
    case 2:
      face[0] = tet[2];
      face[1] = tet[3];
      face[2] = tet[0];
      break;
    case 3:
      face[0] = tet[0];
      face[1] = tet[1];
      face[2] = tet[2];
      break;
    default:
      std::cerr << "ERROR: Invalid face index, code should not reach here!\n";
    }
}

} /* namespace cosmologytools */
