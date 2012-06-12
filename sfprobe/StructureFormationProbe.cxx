#include "StructureFormationProbe.h"

#include "ExtentUtilities.h"
#include "TetrahedronUtilities.h"

// C++ includes
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

    if( this->IsNodeWithinPeriodicBoundaryFringe(tet[node]) )
      {
      return false;
      }

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

  // Mapping of nodes on the langrange mesh to the euler mesh
  std::map<INTEGER,INTEGER> LangrangeID2EulerMeshID;

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
        if(LangrangeID2EulerMeshID.find(lidx)==LangrangeID2EulerMeshID.end())
          {
          this->EulerMesh.Nodes.push_back(nodes[node*3]);
          this->EulerMesh.Nodes.push_back(nodes[node*3+1]);
          this->EulerMesh.Nodes.push_back(nodes[node*3+2]);

          // Euler index
          INTEGER eidx = this->EulerMesh.GetNumberOfNodes()-1;
          LangrangeID2EulerMeshID[ lidx ] = eidx;
          tet[ node ] = eidx;
          } // END if this node has not been mapped
        else
          {
          tet[ node ] = LangrangeID2EulerMeshID[ lidx ];
          } // END else if node has already been mapped

        this->EulerMesh.Connectivity.push_back( tet[ node ] );
        } // END for all tet nodes

      // Compute the volume of the added tetrahedron
      this->Volumes.push_back(
          TetrahedronUtilities::ComputeVolume(
              &nodes[0],&nodes[3],&nodes[6],&nodes[9]) );

      } // END if the tet is mapped succesfully to euler space.

    // The following loop-invariant must hold for correctness
    assert("post: volumes array must be equal to the number of tets!" &&
            this->Volumes.size()==this->EulerMesh.GetNumberOfCells());

    } // END for all langrangian tets

  LangrangeID2EulerMeshID.clear();
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
  // TODO: implement this
}

} /* namespace cosmologytools */
