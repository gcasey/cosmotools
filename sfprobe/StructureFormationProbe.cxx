#include "StructureFormationProbe.h"

#include "ExtentUtilities.h"
#include "TetrahedronUtilities.h"

// C++ includes
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>

namespace cosmologytools {

//------------------------------------------------------------------------------
StructureFormationProbe::StructureFormationProbe()
{
  this->Particles    = NULL;
  this->GlobalIds    = NULL;
  this->Langrange    = NULL;
  this->NumParticles = 0;
  this->Fringe       = 1;
}

//------------------------------------------------------------------------------
StructureFormationProbe::~StructureFormationProbe()
{
  if( this->Langrange != NULL )
    {
    delete this->Langrange;
    }
  this->EulerMesh.Clear();
  this->Volumes.clear();
}

//------------------------------------------------------------------------------
void StructureFormationProbe::SetParticles(
    REAL *particles, INTEGER *globalIds, INTEGER N)
{
  assert("pre: particles array is NULL" && (particles != NULL) );
  assert("pre: global Ids array is NULL" && (globalIds != NULL) );
  assert("pre: Number of particles must be greater than 1" && (N >= 1) );

  this->EulerMesh.Clear();

  this->NumParticles = N;
  this->Particles = particles;
  this->GlobalIds = globalIds;

  for(INTEGER pos=0; pos < N; ++pos )
    {
    this->Global2PositionMap[ globalIds[pos] ] = pos;
    } // END for all particle positions
}


//------------------------------------------------------------------------------
void StructureFormationProbe::BuildLangrangianMesh(
        REAL origin[3], REAL spacing[3], INTEGER extent[6])
{
  if( this->Langrange != NULL )
    {
    delete this->Langrange;
    }

  this->Langrange = new LangrangianTesselator();
  this->Langrange->SetOrigin(origin);
  this->Langrange->SetSpacing(spacing);
  this->Langrange->SetGridExtent(extent);
  this->Langrange->Tesselate();
}

//------------------------------------------------------------------------------
LangrangianTesselator* StructureFormationProbe::GetLangrangeTesselator()
{
  return this->Langrange;
}

//------------------------------------------------------------------------------
bool StructureFormationProbe::IsNodeWithinFringeBounds(INTEGER nodeIdx)
{
  assert("pre: Langrange tesselator object is NULL" &&
         (this->Langrange != NULL) );

  INTEGER idx = this->Global2PositionMap[ nodeIdx ];

  REAL pnt[3];
  pnt[0] = this->Particles[idx*3];
  pnt[1] = this->Particles[idx*3+1];
  pnt[2] = this->Particles[idx*3+2];

  std::cout << "FRINGE: " << this->Fringe << std::endl;
  std::cout.flush();

  REAL iBounds[6];
  this->Langrange->GetBounds(iBounds,this->Fringe);

  return(
      (pnt[0] >= iBounds[0]) && (pnt[0] <= iBounds[1]) &&
      (pnt[1] >= iBounds[2]) && (pnt[1] <= iBounds[3]) &&
      (pnt[2] >= iBounds[4]) && (pnt[2] <= iBounds[5])
      );
}

//------------------------------------------------------------------------------
bool StructureFormationProbe::IsNodeWithinPeriodicBoundaryFringe(
          INTEGER nodeIdx )
{
  assert("pre: Langrange tesselator object is NULL" &&
          (this->Langrange != NULL) );

  INTEGER idx = this->Global2PositionMap[ nodeIdx ];

  REAL pnt[3];
  pnt[0] = this->Particles[idx*3];
  pnt[1] = this->Particles[idx*3+1];
  pnt[2] = this->Particles[idx*3+2];

  REAL bounds[6];
  this->Langrange->GetBounds(bounds);
  for( int i=0; i < 3; ++i )
    {
    REAL L     = bounds[i*2+1]-bounds[i*2];
    REAL delta = static_cast<REAL>( (L*this->Fringe) )/100.0;
    if( pnt[i] < bounds[i*2]+delta ||
        pnt[i] > bounds[i*2+1]-delta)
      {
      return true;
      } // END if
    } // END for
  return false;
}

//------------------------------------------------------------------------------
bool StructureFormationProbe::MapTetToEulerSpace(
        INTEGER tet[4], REAL nodes[12])
{
  for(int node=0; node < 4; ++node )
    {
    if( this->Global2PositionMap.find(tet[node])==
            this->Global2PositionMap.end() )
      {
      return false;
      }

    if( !this->IsNodeWithinFringeBounds(tet[node]))
      {
      return false;
      }
//    if( this->IsNodeWithinPeriodicBoundaryFringe(tet[node]) )
//      {
//      return false;
//      }

    INTEGER ppos = this->Global2PositionMap[ tet[node] ];
    for( int dim=0; dim < 3; ++dim )
      {
      nodes[ node*3+dim ] = this->Particles[ ppos*3+dim ];
      }

    } // END for all nodes
  return true;
}

//------------------------------------------------------------------------------
void StructureFormationProbe::BuildEulerMesh()
{
  // Sanity check
  assert(
    "pre: Must construct a langrangian mesh first before an euler mesh!" &&
    (this->Langrange != NULL) );

  this->LangrangeNode2EulerNode.clear();
  this->LangrangeTet2EulerTet.clear();
  this->EulerMesh.Clear();
  this->Volumes.clear();

  // Reserve space to store euler mesh
  this->EulerMesh.Stride = 4;

  // Reserve space for Euler mesh
  this->EulerMesh.Nodes.reserve(
      ExtentUtilities::ComputeNumberOfNodes(
            const_cast<INTEGER*>(this->Langrange->GetExtent())));
  this->EulerMesh.Connectivity.reserve(
      this->EulerMesh.Stride*this->Langrange->GetNumTets());
  this->Volumes.reserve(this->Langrange->GetNumTets());

  for(INTEGER idx=0; idx < this->Langrange->GetNumTets(); ++idx)
    {
    INTEGER tet[4];
    this->Langrange->GetTetConnectivity(idx,tet);

    REAL nodes[12];
    if( this->MapTetToEulerSpace(tet,nodes) )
      {
      // Insert vertices & tetrahedra in the mesh
      for( int node=0; node < 4; ++node )
        {
        // Langrange index
        INTEGER lidx = tet[ node ];
        if(this->LangrangeNode2EulerNode.find(lidx)==
           this->LangrangeNode2EulerNode.end())
          {
          this->EulerMesh.Nodes.push_back(nodes[node*3]);
          this->EulerMesh.Nodes.push_back(nodes[node*3+1]);
          this->EulerMesh.Nodes.push_back(nodes[node*3+2]);

          // Euler index
          INTEGER eidx = this->EulerMesh.GetNumberOfNodes()-1;
          this->LangrangeNode2EulerNode[ lidx ] = eidx;
          tet[ node ] = eidx;
          } // END if this node has not been mapped
        else
          {
          tet[ node ] = this->LangrangeNode2EulerNode[ lidx ];
          } // END else if node has already been mapped

        this->EulerMesh.Connectivity.push_back( tet[ node ] );
        } // END for all tet nodes

      this->LangrangeTet2EulerTet[idx]=this->EulerMesh.GetNumberOfCells()-1;

      // Compute the volume of the added tetrahedron
      this->Volumes.push_back(
          TetrahedronUtilities::ComputeVolume(
              &nodes[0],&nodes[3],&nodes[6],&nodes[9]) );

      } // END if the tet is mapped succesfully to euler space.

    // The following loop-invariant must hold for correctness
    assert("post: volumes array must be equal to the number of tets!" &&
            this->Volumes.size()==this->EulerMesh.GetNumberOfCells());

    } // END for all langrangian tets

}

//------------------------------------------------------------------------------
INTEGER StructureFormationProbe::GetNumberOfStreams(REAL pnt[3])
{
  // TODO: implement this
  return 0;
}

//------------------------------------------------------------------------------
void StructureFormationProbe::GetEulerMesh(
        std::vector<REAL> &nodes, std::vector<INTEGER> &tets,
        std::vector<REAL> &volumes)
{
  // Copy internal data-structures to user-supplied vectors
  nodes   = this->EulerMesh.Nodes;
  tets    = this->EulerMesh.Connectivity;
  volumes = this->Volumes;
}

//------------------------------------------------------------------------------
void StructureFormationProbe::ExtractCausticSurfaces(
        std::vector<REAL> &nodes, std::vector<INTEGER> &triangles)
{
  // Sanity check
  assert("pre: Must construct a langrangian mesh first before an euler mesh!"
         && (this->Langrange != NULL) );

  // STEP 0: Initialize supplied vectors
  nodes.clear();
  triangles.clear();

  // STEP 1: Get Langrange Mesh faces
  std::vector< INTEGER > faces;
  this->Langrange->GetFaces( faces );

  // STEP 2: Loop through Langrange mesh faces and determine if they correspond
  // to a face on the caustics surface mesh.
  std::vector<INTEGER> tets;
  INTEGER face[3];
  for(INTEGER fidx=0; fidx < faces.size()/3; ++fidx)
    {
    for( int node=0; node < 3; ++node )
      {
      face[node] = faces[fidx*3+node];
      } // END for all face nodes

    this->Langrange->GetAdjacentTets(face,tets);
    assert("pre: face has more than 2 adjacent tets!" &&
            (tets.size() >= 1) && (tets.size() <= 2) );

    if( tets.size() != 2 )
      {
      // This is a boundary face, skip!
      continue;
      }

    if( (this->LangrangeTet2EulerTet.find(tets[0])==
         this->LangrangeTet2EulerTet.end()) ||
        (this->LangrangeTet2EulerTet.find(tets[1]) ==
         this->LangrangeTet2EulerTet.end()) )
      {
      // This face abutts a tet that wasn't mapped to the euler mesh, skip!
      continue;
      }

    // Get the tet indices in the euler mesh
    INTEGER tetIdx1 = this->LangrangeTet2EulerTet[tets[0]];
    INTEGER tetIdx2 = this->LangrangeTet2EulerTet[tets[1]];

    std::map<INTEGER,INTEGER> eulerMesh2CausticSurface;

    if( this->VolumesHaveOppositeSigns(tetIdx1,tetIdx2) )
      {
      INTEGER t[3];
      for( int tnode=0; tnode < 3; ++tnode)
        {
        assert("caustics surface mesh node cannot be found in Euler mesh!" &&
                this->LangrangeNode2EulerNode.find(face[tnode]) !=
                this->LangrangeNode2EulerNode.end());

        // Get the node index w.r.t. the euler mesh
       INTEGER nodeIdx = this->LangrangeNode2EulerNode[face[tnode]];

        if( eulerMesh2CausticSurface.find(nodeIdx) ==
             eulerMesh2CausticSurface.end() )
          {
          nodes.push_back(this->EulerMesh.Nodes[nodeIdx*3 ]);
          nodes.push_back(this->EulerMesh.Nodes[nodeIdx*3+1]);
          nodes.push_back(this->EulerMesh.Nodes[nodeIdx*3+2]);
          eulerMesh2CausticSurface[nodeIdx] = (nodes.size()/3)-1;
          } // END if
        t[tnode] = eulerMesh2CausticSurface[nodeIdx];
        triangles.push_back( t[tnode] );
        } // END for all tnodes

      } // END if the volumes of the adjacent tets is opposite

    } // END for all mesh faces
}

//------------------------------------------------------------------------------
bool StructureFormationProbe::VolumesHaveOppositeSigns(
        INTEGER tetIdx1, INTEGER tetIdx2)
{
  assert("pre: tet index 1 is out-of-bounds!" &&
     (tetIdx1 >= 0) && (tetIdx1 < this->Volumes.size() ) );
  assert("pre: tet index 2 is out-of-bounds!" &&
     (tetIdx2 >= 0) && (tetIdx2 < this->Volumes.size() ) );

  REAL v1 = this->Volumes[tetIdx1];
  REAL v2 = this->Volumes[tetIdx2];
  return( ((v1*v2) < 0.0f) );
}

} /* namespace cosmologytools */
