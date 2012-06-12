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
