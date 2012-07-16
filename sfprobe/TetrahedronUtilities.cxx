#include "TetrahedronUtilities.h"

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

} /* namespace cosmologytools */
