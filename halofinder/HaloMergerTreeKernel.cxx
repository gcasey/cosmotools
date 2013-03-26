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
  this->MergerTreeThreshold = 50;
  this->ZombieCutOff = 5;
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

  this->DeadHalos.clear();
  this->SplitHalos.clear();
  this->MergeHalos.clear();

  int nrows = this->Sizes[0];
  int ncol  = this->Sizes[1];
  this->HaloSimilarityMatrix.resize( nrows*ncol );

  this->MatrixColumnSum.resize(ncol,0);

  int rowCount = 0;
  for( int row=0; row < nrows; ++row )
    {
    // STEP 0: Initialize rowCount to 0
    rowCount = 0;

    // STEP 1: Propagate zombies across timesteps without further checking
    if( this->Halos1[row].HaloType == ZOMBIEHALO &&
        this->Halos1[row].Count > this->ZombieCutOff )
      {
      this->DeadHalos.insert(row);
      continue;
      }

    // STEP 2: Loop through all columns in this row
    for( int col=0; col < ncol; ++col )
      {
      int overlap = this->Halos1[row].Intersect(&this->Halos2[col]);
      this->HaloSimilarityMatrix[row*ncol+col] = overlap;
      if( overlap >= this->MergerTreeThreshold )
        {
        this->MatrixColumnSum[ col ]++;
        ++rowCount;
        } // END if
      } // END for all columns

     // STEP 3: Detect 'DEATH' or 'SPLIT" event
     switch( rowCount )
       {
       case 0:
         // DEATH event
         this->DeadHalos.insert( row );
         break;
       case 1:
         // CONTINUATION or MERGE event. The column sum matrix
         // must be examined to determine which event it is.
         break;
       default:
         this->SplitHalos.insert( row );
       } // END switch
    } // END for all rows

  // Loop through the matrix-column sum & detect merge
  unsigned int colIdx=0;
  for(;colIdx < this->MatrixColumnSum.size(); colIdx++)
    {
    switch( this->MatrixColumnSum[colIdx] )
      {
      case 0:
        // BIRTH event
        break;
      case 1:
        // CONTINUATION event
        break;
      default:
        // MERGE event
        this->MergeHalos.insert( colIdx );
      } // END switch
    } // END for all columns
}

//------------------------------------------------------------------------------
void HaloMergerTreeKernel::UpdateHaloEvolutionTree(
        DistributedHaloEvolutionTree *t)
{
  assert("pre: halo evolution tree should not be NULL" && (t!=NULL) );

  // STEP 1: Create Edges
  int nrows = this->Sizes[0];
  int ncol  = this->Sizes[1];
  for(int row=0; row < nrows; ++row)
    {

    if( this->IsHaloDead( row )  )
      {
      continue;
      }

    for(int col=0; col < ncol; ++col)
      {

      int overlap = this->HaloSimilarityMatrix[row*ncol+col];
      if(overlap >= this->MergerTreeThreshold)
        {
        int event = MergerTreeEvent::UNDEFINED;
        if( this->IsHaloMerge(col) )
          {
          if( this->Halos1[row].HaloType == ZOMBIEHALO )
            {
            event = MergerTreeEvent::REBIRTH;
            }
          else
            {
            event = MergerTreeEvent::MERGE;
            }
          }
        else if(this->IsHaloSplit(row))
          {
          event = MergerTreeEvent::SPLIT;
          }
        else
          {
          if( this->Halos1[row].HaloType == ZOMBIEHALO )
            {
            event = MergerTreeEvent::REBIRTH;
            this->Halos1[row].Count = 0;
            }
          else
            {
            event = MergerTreeEvent::CONTINUATION;
            }
          }

        t->InsertNode(this->Halos1[row]);
        t->InsertNode(this->Halos2[col]);
        t->CreateEdge(
            this->Halos1[row].GetHashCode(),
            this->Halos2[col].GetHashCode(),
            overlap,
            event);
        } // if the halos are similar
      else
        {
//        std::cout << "NOTE: Deal with nodes that fall below the threshold!\n";
//        std::cout.flush();
        }

      } // END for all columns

    } // END for all rows

}

} /* namespace cosmotk */
