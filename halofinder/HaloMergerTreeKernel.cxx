#include "HaloMergerTreeKernel.h"

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
  this->Verbose = true;
  this->NumberOfBirths = 0;
  this->NumberOfRebirths = 0;
}

//------------------------------------------------------------------------------
HaloMergerTreeKernel::~HaloMergerTreeKernel()
{
  this->HaloSimilarityMatrix.clear();
  this->MatrixColumnSum.clear();
  this->MatrixRowSum.clear();
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

  // STEP 0: Initialize data-structures
  this->DeadHalos.clear();
  this->SplitHalos.clear();
  this->MergeHalos.clear();
  this->ProcessedColumns.clear();
  this->NumberOfBirths = 0;
  this->NumberOfRebirths = 0;

  int nrows = this->Sizes[0];
  int ncol  = this->Sizes[1];
  this->HaloSimilarityMatrix.resize( nrows*ncol );
  this->MatrixRowSum.resize(nrows,0);
  this->MatrixColumnSum.resize(ncol,0);

  // STEP 1: Compute similarity matrix
  for( int row=0; row < nrows; ++row )
    {
    // STEP 1.0: Propagate zombies across timesteps without further checking
    // if the given cut-off criteria is met, i.e., if a zombie is persistently
    // a zombie for more that ZombieCutOff time-steps then, it is a zombie and
    // we don't bother checking it any more.
    if( this->Halos1[row].HaloType == ZOMBIEHALO &&
        this->Halos1[row].Count > this->ZombieCutOff )
      {
      this->DeadHalos.insert(row);
      continue;
      }

    // STEP 1.1: Loop through all columns in this row
    for( int col=0; col < ncol; ++col )
      {
      int overlap = this->Halos1[row].Intersect(&this->Halos2[col]);
      this->HaloSimilarityMatrix[row*ncol+col] = overlap;
      if( this->MajorityRuleCheck(overlap) )
        {
        this->MatrixColumnSum[ col ]++;
        this->MatrixRowSum[ row ]++;
        } // END if
      } // END for all columns
    } // END for all rows
}

//------------------------------------------------------------------------------
void HaloMergerTreeKernel::UpdateHaloEvolutionTree(
        DistributedHaloEvolutionTree *mergerTree)
{
  assert("pre: halo evolution tree should not be NULL" &&
         (mergerTree !=NULL) );

  int nrows    = this->Sizes[0];
  int rowEvent = MergerTreeEvent::UNDEFINED;
  for(int row=0; row < nrows; ++row)
    {
    // STEP 0: If this halo has been permanently designated as a zombie, skip it.
    if( this->IsHaloDead(row) )
      {
      continue;
      }
    rowEvent = MergerTreeEvent::UNDEFINED;

    Halo* prevHalo = &this->Halos1[row];
    assert("pre: previous halo is NULL!" && (prevHalo != NULL) );

    // STEP 1: If the halo is not already in the halo evolution tree, insert it
    // and mark it as a birth. For example, this happens when processing the
    // first time-step pair.
    if(!mergerTree->HasNode(prevHalo->GetHashCode()))
      {
      rowEvent = MergerTreeEvent::BIRTH;
      ++this->NumberOfBirths;
      unsigned char bitmask;
      MergerTreeEvent::Reset(bitmask);
      MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::BIRTH);
      this->InsertHalo(prevHalo,bitmask,mergerTree);
      }

    // STEP 2: Check if the halo in the previous time-step died in the next
    // time-step.
    if(this->MatrixRowSum[row] == 0)
      {
      rowEvent = MergerTreeEvent::DEATH;
      this->DeadHalos.insert(row);
      continue;
      }
    // STEP 3: Check if this is a continuation or rebirth
    else if( this->MatrixRowSum[row] == 1)
      {
      // Do nothing here, fall through to detect event in the next time-step
      rowEvent = MergerTreeEvent::CONTINUATION;
      }
    // STEP 3: Check if the halo is split in the next time-step
    else if(this->MatrixRowSum[row]  > 1)
      {
      rowEvent = MergerTreeEvent::SPLIT;
      this->SplitHalos.insert(row);
      }
    else
      {
      assert("ERROR: Code should not reach here!" && false);
      }

    // STEP 4: Detect event for the halo of the previous time-step, to the
    // current timestep.
    this->DetectEvent(row,rowEvent,prevHalo,mergerTree);
    } // END for all rows

  if( this->Verbose )
    {
    std::cout << "NUMBER OF DEATHS: " << this->DeadHalos.size()  << std::endl;
    std::cout << "NUMBER OF SPLITS: " << this->SplitHalos.size() << std::endl;
    std::cout << "NUMBER OF MERGES: " << this->MergeHalos.size() << std::endl;
    std::cout << "NUMBER OF BIRTHS: " << this->NumberOfBirths << std::endl;
    std::cout << "=========\n";
    std::cout.flush();
    }
}

//------------------------------------------------------------------------------
void HaloMergerTreeKernel::DetectEvent(
        const int row, const int rowEvent,
        Halo* prevHalo, DistributedHaloEvolutionTree *mergerTree)
{
  assert("pre: previous halo pointer is NULL!" &&
         (prevHalo != NULL) );
  assert("pre: merger-tree pointer is NULL!" &&
         (mergerTree != NULL) );
  assert("pre: rowEvent is undefined!" &&
         (rowEvent != MergerTreeEvent::UNDEFINED) );

  int ncol     = this->Sizes[1];
  int colEvent = MergerTreeEvent::UNDEFINED;
  for(int col=0; col < ncol; ++col)
    {
    colEvent = MergerTreeEvent::UNDEFINED;

    // STEP 0: Short-circuit if we have finished processing this column
    if( this->ProcessedColumns.find(col) != this->ProcessedColumns.end() )
      {
      continue;
      }

    // STEP 1: Get the current halo and percent overlap computed in
    // ComputeMergerTree
    Halo* currHalo = &this->Halos2[col];
    assert("pre: current halo is NULL!" && (currHalo != NULL) );
    assert("pre: halo at current timestep cannot be a zombie" &&
           (currHalo->HaloType != ZOMBIEHALO));
    int overlap = this->HaloSimilarityMatrix[row*ncol+col];

    // STEP 2: Check if this is a new halo
    if( this->MatrixColumnSum[col] == 0)
      {
      // This is a birth of a new halo that is not related to any halos from
      // the previous time-step
      ++this->NumberOfBirths;
      unsigned char bitmask;
      MergerTreeEvent::Reset(bitmask);
      MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::BIRTH);
      this->InsertHalo(currHalo,bitmask,mergerTree);
      ProcessedColumns.insert(col);
      }
    // STEP 3: Check if the majority rule passes. If it passes, detect the
    // event type and insert the node to the tree, linking it to the halo
    // from the previous time-step.
    else if(this->MajorityRuleCheck(overlap))
      {
      assert("pre: previous halo not in tree!" &&
              mergerTree->HasNode(prevHalo->GetHashCode()));

      unsigned char bitmask;
      MergerTreeEvent::Reset(bitmask);

      switch(this->MatrixColumnSum[col])
        {
        case 1:
          {
          MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::CONTINUATION);

          if( prevHalo->HaloType==ZOMBIEHALO )
            {
            MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::REBIRTH);
            ++this->NumberOfRebirths;
            }

          if( rowEvent == MergerTreeEvent::SPLIT )
            {
            MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::SPLIT);
            }

          this->InsertHalo(currHalo,bitmask,mergerTree);
          mergerTree->LinkHalos(
              prevHalo->GetHashCode(),currHalo->GetHashCode());
          this->ProcessedColumns.insert( col );
          }
          break;
        default:
          assert("pre: merger must be with 2 or more halos" &&
                  (this->MatrixColumnSum[col] >= 2) );

          this->MergeHalos.insert(col);

          MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::MERGE);

          if( prevHalo->HaloType==ZOMBIEHALO )
            {
            MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::REBIRTH);
            }

          if( rowEvent == MergerTreeEvent::SPLIT )
            {
            MergerTreeEvent::SetEvent(bitmask,MergerTreeEvent::SPLIT);
            }

          this->InsertHalo(currHalo,bitmask,mergerTree);
          mergerTree->LinkHalos(
              prevHalo->GetHashCode(),currHalo->GetHashCode());
        } // END switch
      } // END if majority rule passes
    } // END for all columns

}

} /* namespace cosmotk */
