#include "HaloMergerTreeKernel.h"
#include "Halo.h"

#include <cassert>

namespace cosmotk
{

HaloMergerTreeKernel::HaloMergerTreeKernel()
{
  this->MergerTreeThreshold = 10;
}

//------------------------------------------------------------------------------
HaloMergerTreeKernel::~HaloMergerTreeKernel()
{
  this->HaloSimilarityMatrix.clear();
}


//------------------------------------------------------------------------------
void HaloMergerTreeKernel::ComputeMergerTree(
    const int t1, std::vector< Halo > &haloSet1,
    const int t2, std::vector< Halo > &haloSet2 )
{
  assert("pre: t2 > t1 must hold true" && (t2 > t1) );

  int nrows = static_cast<int>( haloSet1.size() );
  int ncol  = static_cast<int>( haloSet2.size() );
  this->HaloSimilarityMatrix.resize( nrows*ncol );

  for( int row=0; row < nrows; ++row )
    {
    for( int col=0; col < ncol; ++col )
      {
      int overlap = haloSet1[ row ].Intersect( &haloSet2[col ] );
      this->HaloSimilarityMatrix[row*ncol+col] = overlap;
      } // END for all columns
    } // END for all rows
}

//------------------------------------------------------------------------------
void HaloMergerTreeKernel::UpdateHaloEvolutionTree()
{
  // TODO: implement this

}

} /* namespace cosmotk */
