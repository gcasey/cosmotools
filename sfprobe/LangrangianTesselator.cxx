#include "LangrangianTesselator.h"

#include "ExtentUtilities.h"

// C/C++ includes
#include <cstdlib> // For definition of NULL
#include <cstring> // For memcpy
#include <cassert> // For assert

namespace cosmologytools
{

LangrangianTesselator::LangrangianTesselator()
{
  this->Initialize();
}

//-----------------------------------------------------------------------------
LangrangianTesselator::~LangrangianTesselator()
{
  if( this->Connectivity != NULL )
    {
    delete [] this->Connectivity;
    }
  this->Initialize();
}

//-----------------------------------------------------------------------------
void LangrangianTesselator::SetGridExtent(int extent[6])
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
void LangrangianTesselator::Initialize()
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
void LangrangianTesselator::GetPoint(INTEGER ijk[3], REAL pnt[3])
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
 ExtentUtilities::GetStructuredCoordinates(nodeIdx,ijk,this->Extent); \
 this->GetPoint( ijk, V);                                              \
}

//------------------------------------------------------------------------------
void LangrangianTesselator::GetLangrangianTet(
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
void LangrangianTesselator::GetTetConnectivity(
        const INTEGER tetIdx, INTEGER tet[4])
{
  assert("pre: tetIdx is out-of-bounds" && (tetIdx < this->NumTets) );

  for( int i=0; i < 4; ++i )
    {
    tet[i] = this->Connectivity[tetIdx*4+i];
    }
}

//------------------------------------------------------------------------------
#define ADDTET(tetIdx, V0, V1, V2, V3 )   \
    this->Connectivity[tetIdx*4]   = V0;  \
    this->Connectivity[tetIdx*4+1 ]= V1;  \
    this->Connectivity[tetIdx*4+2 ]= V2;  \
    this->Connectivity[tetIdx*4+3 ]= V3;  \
    ++tetIdx;

//-----------------------------------------------------------------------------
void LangrangianTesselator::BuildTesselation()
{
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
          ADDTET(tetIdx,V[1],V[2],V[3],V[6]);
          ADDTET(tetIdx,V[0],V[1],V[3],V[4]);
          ADDTET(tetIdx,V[4],V[6],V[5],V[1]);
          ADDTET(tetIdx,V[4],V[6],V[7],V[3]);
          ADDTET(tetIdx,V[3],V[6],V[4],V[1]);
          }
        else
          {
          ADDTET(tetIdx,V[1],V[5],V[2],V[0]);
          ADDTET(tetIdx,V[2],V[3],V[0],V[7]);
          ADDTET(tetIdx,V[2],V[5],V[6],V[7]);
          ADDTET(tetIdx,V[0],V[7],V[4],V[5]);
          ADDTET(tetIdx,V[0],V[2],V[7],V[5]);
          }

        } // END for all k
      } // END for all j
    } // END for all i
}

} /* namespace cosmologytools */
