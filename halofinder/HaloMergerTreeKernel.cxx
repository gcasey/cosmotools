#include "HaloMergerTreeKernel.h"

#include "DistributedHaloEvolutionTree.h"
#include "Halo.h"

#include <cassert>

namespace cosmotk
{

HaloMergerTreeKernel::HaloMergerTreeKernel()
{
  this->Timesteps[0] = this->Timesteps[1] = 0;
  this->Sizes[0]     = this->Sizes[1]     = 1;
  this->Halos1 = NULL;
  this->Halos2 = NULL;
  this->MergerTreeThreshold = 10;
}

//------------------------------------------------------------------------------
HaloMergerTreeKernel::~HaloMergerTreeKernel()
{
  this->HaloSimilarityMatrix.clear();
}


//------------------------------------------------------------------------------
void HaloMergerTreeKernel::UpdateMergerTree(
      const int t1, Halo *haloSet1, const int M,
      const int t2, Halo *haloSet2, const int N,
      DistributedHaloEvolutionTree *t )
{
  this->RegisterHalos( t1,haloSet1,M, t2,haloSet2,N );
  this->ComputeMergerTree();
  this->UpdateHaloEvolutionTree( t );
}

//------------------------------------------------------------------------------
void HaloMergerTreeKernel::RegisterHalos(
    const int t1, Halo *haloSet1, const int M,
    const int t2, Halo *haloSet2, const int N)
{
  assert("pre: t1 < t2" && (t1 < t2) );
  assert("pre: halo set 1 should not be empty" && (M >= 1) );
  assert("pre: halo set 2 should not be empty" && (N >= 1) );
  assert("pre: null haloset!" && (haloSet1 != NULL) );
  assert("pre: null haloset!" && (haloSet2 != NULL) );

  this->Timesteps[0] = t1;
  this->Timesteps[1] = t2;
  this->Sizes[0] = M;
  this->Sizes[1] = N;
  this->Halos1 = haloSet1;
  this->Halos2 = haloSet2;
}

//------------------------------------------------------------------------------
void HaloMergerTreeKernel::ComputeMergerTree( )
{
  assert("pre: halos != NULL" && (this->Halos1 != NULL) );
  assert("pre: halos != NULL" && (this->Halos2 != NULL) );

  int nrows = this->Sizes[0];
  int ncol  = this->Sizes[1];
  this->HaloSimilarityMatrix.resize( nrows*ncol );

  for( int row=0; row < nrows; ++row )
    {
    for( int col=0; col < ncol; ++col )
      {
      int overlap = this->Halos1[row].Intersect(&this->Halos2[col]);
      this->HaloSimilarityMatrix[row*ncol+col] = overlap;
      } // END for all columns
    } // END for all rows
}

//------------------------------------------------------------------------------
void HaloMergerTreeKernel::UpdateHaloEvolutionTree(
        DistributedHaloEvolutionTree *t)
{
  assert("pre: halo evolution tree should not be NULL" && (t!=NULL) );

  // STEP 0: Append nodes to the tree
  if( t->IsEmpty() )
    {
    t->AppendNodes(this->Halos1,this->Sizes[0]);
    t->AppendNodes(this->Halos2,this->Sizes[1]);
    }
  else
    {
    t->AppendNodes(this->Halos2,this->Sizes[1]);
    }

  // STEP 1: Create Edges
  int nrows = this->Sizes[0];
  int ncol  = this->Sizes[1];
  for(int row=0; row < nrows; ++row)
    {
    for(int col=0; col < ncol; ++col)
      {
      int overlap = this->HaloSimilarityMatrix[row*ncol+col];
      if(overlap > this->MergerTreeThreshold)
        {
        t->CreateEdge(
            this->Halos1[row].GetHashCode(),
            this->Halos2[col].GetHashCode(),
            overlap);
        } // if the halos are similar
      } // END for all columns
    } // END for all rows

}

} /* namespace cosmotk */
