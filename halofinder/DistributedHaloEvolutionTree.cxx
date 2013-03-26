#include "DistributedHaloEvolutionTree.h"

// CosmologyTools includes
#include "Halo.h"

// STL includes
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <set>
#include <limits>

//------------------------------------------------------------------------------

namespace cosmotk
{

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::DistributedHaloEvolutionTree()
{
  this->Communicator = MPI_COMM_NULL;
  this->FileName = "HaloEvolutionTree.dat";
  this->IOFormat = MergerTreeFileFormat::GENERIC_IO;
}

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::~DistributedHaloEvolutionTree()
{
  this->Clear();
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::InsertNode(Halo &halo)
{
  this->Nodes[ halo.GetHashCode() ] = halo;
  this->Nodes[ halo.GetHashCode() ].ParticleIds.clear();
  if( halo.HaloType != GHOSTHALO )
    {
    this->UpdateNodeCounter( halo.TimeStep );
    }

  if( halo.HaloType == ZOMBIEHALO )
    {
    this->UpdateZombieCounter( halo.TimeStep );
    }
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::AppendNodes(
      Halo *halos, const int N)
{
  for( int halo=0; halo < N; ++halo )
    {
    this->Nodes[ halos[ halo ].GetHashCode() ] = halos[ halo ];

    if( halos[ halo ].HaloType != GHOSTHALO )
      {
      this->UpdateNodeCounter( halos[ halo ].TimeStep );
      }

    if( halos[ halo ].HaloType == ZOMBIEHALO )
      {
      this->UpdateZombieCounter( halos[ halo ].TimeStep );
      }

    // Remove the particle Ids of the object
    this->Nodes[ halos[ halo ].GetHashCode() ].ParticleIds.clear();
    } // END for all halos
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::CreateEdge(
      const std::string halo1,
      const std::string halo2,
      int w,
      int e)
{
  assert("pre: node halo1 does not exist" && this->HasNode(halo1) );
  assert("pre: node halo2 does not exist" && this->HasNode(halo2) );
  this->Edges.push_back( halo1 );
  this->Edges.push_back( halo2 );
  this->EdgeWeights.push_back( w );
  this->EdgeEvents.push_back( e );

  if( this->NodeDescendants.find(halo1) != this->NodeDescendants.end() )
    {
    this->NodeDescendants[ halo1 ].insert( halo2 );
    }
  else
    {
    std::set< std::string > descendants;
    descendants.insert( halo2 );
    this->NodeDescendants[ halo1 ] = descendants;
    }
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::Clear()
{
  this->Nodes.clear();
  this->Node2UniqueIdx.clear();
  this->Edges.clear();
  this->EdgeWeights.clear();
  this->EdgeEvents.clear();
}

//------------------------------------------------------------------------------
bool DistributedHaloEvolutionTree::IsEmpty()
{
  return( this->Nodes.empty() );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfNodes(const int tstep)
{
  int N = 0;
  if(this->NodeCounter.find(tstep) != this->NodeCounter.end() )
    {
    N = this->NodeCounter[ tstep ];
    }
  else
    {
    std::cerr << "WARNING: No tree nodes for requrested time-step!\n";
    }
  return( N );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfZombieNodes(
    const int tstep)
{
  int N = 0;
  if( this->ZombieCounter.find(tstep) != this->ZombieCounter.end() )
    {
    N = this->ZombieCounter[ tstep ];
    }
  return( N );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfNodes()
{
  return( this->Nodes.size() );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfEdges()
{
  return( this->EdgeWeights.size() );
}

//------------------------------------------------------------------------------
bool DistributedHaloEvolutionTree::HasNode(std::string hashCode)
{
  if( this->Nodes.find(hashCode) == this->Nodes.end() )
    {
    return false;
    }
  return true;
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::UpdateNodeCounter(const int tstep)
{
  if( this->NodeCounter.find(tstep) != this->NodeCounter.end() )
    {
    this->NodeCounter[tstep]++;
    }
  else
    {
    this->NodeCounter[tstep] = 1;
    }
}

//-----------------------------------------------------------------------------
void DistributedHaloEvolutionTree::UpdateZombieCounter(const int tstep)
{
  if( this->ZombieCounter.find(tstep) != this->ZombieCounter.end() )
    {
    this->ZombieCounter[tstep]++;
    }
  else
    {
    this->ZombieCounter[tstep] = 1;
    }
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::GetNodeRangeForProcess(ID_T r[2])
{
  assert("pre: NULL communicator!" && (this->Communicator != MPI_COMM_NULL));

  // STEP 0: Get number of nodes in this process
  ID_T numNodes = static_cast<ID_T>( this->Nodes.size() );

  // STEP 1: Do a prefix sum
  ID_T mySum = -1;
  MPI_Scan(&numNodes,&mySum,1,MPI_ID_T,MPI_SUM,this->Communicator);

  // STEP 2: Get range in this process
  r[0] = (mySum-1)-numNodes+1;
  r[1] = mySum-1;
}

//------------------------------------------------------------------------------
ID_T DistributedHaloEvolutionTree::GetNodeIndex(std::string hashCode)
{
  assert( "pre: Cannot find Unique index for node!" &&
      this->Node2UniqueIdx.find(hashCode)!=this->Node2UniqueIdx.end());
  return(this->Node2UniqueIdx[hashCode]);
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::RelabelTreeNodes()
{
  assert("pre: NULL communicator!" && (this->Communicator != MPI_COMM_NULL));

  // STEP 0: Get node range on this process
  ID_T range[2];
  this->GetNodeRangeForProcess(range);
  assert( "pre: Assigned range does not match number of nodes!" &&
      ( (range[1]-range[0]+1)==this->GetNumberOfNodes() ));

  // STEP 1: Assign unique indices to nodes
  std::map<std::string,Halo>::iterator iter = this->Nodes.begin();
  for(ID_T idx=range[0]; idx <= range[1]; ++idx, ++iter)
    {
    this->Node2UniqueIdx[ iter->first ] = idx;
    } // END for the entire range

  // STEP 2: Resolve duplicates ?
  // TODO: implement this
}

//------------------------------------------------------------------------------
std::string DistributedHaloEvolutionTree::ToString()
{
  this->RelabelTreeNodes();

  std::ostringstream oss;
  oss << "NUMNODES " << this->GetNumberOfNodes() << std::endl;
  std::map<std::string, Halo>::iterator NodeIter = this->Nodes.begin();
  oss << "ID\t";
  oss << "TIMESTEP\t";
  oss << "TAG\t";
  oss << "REDSHIFT\t";
  oss << "CENTER_X\t";
  oss << "CENTER_Y\t";
  oss << "CENTER_Z\t";
  oss << "MEAN_CENTER_X\t";
  oss << "MEAN_CENTER_Y\t";
  oss << "MEAN_CENTER_Z\t";
  oss << "V_X\t";
  oss << "V_Y\t";
  oss << "V_Z\t";
  oss << "MASS\t";
  oss << "DESCENDANTS";
  oss << std::endl;
  for(;NodeIter != this->Nodes.end(); ++NodeIter)
    {
    oss << std::scientific
        << std::setprecision(std::numeric_limits<POSVEL_T>::digits10);
    oss << this->GetNodeIndex(NodeIter->first) << "\t";
    oss << NodeIter->second.TimeStep  << "\t";
    oss << NodeIter->second.Tag       << "\t";
    oss << NodeIter->second.Redshift  << "\t";
    oss << NodeIter->second.Center[0] << "\t";
    oss << NodeIter->second.Center[1] << "\t";
    oss << NodeIter->second.Center[2] << "\t";
    oss << NodeIter->second.MeanCenter[0] << "\t";
    oss << NodeIter->second.MeanCenter[1] << "\t";
    oss << NodeIter->second.MeanCenter[2] << "\t";
    oss << NodeIter->second.AverageVelocity[0] << "\t";
    oss << NodeIter->second.AverageVelocity[1] << "\t";
    oss << NodeIter->second.AverageVelocity[2] << "\t";
    oss << NodeIter->second.HaloMass << "\t";

    if( this->NodeDescendants.find(NodeIter->first) !=
        this->NodeDescendants.end())
      {
      std::set< std::string >::iterator descentIter =
          this->NodeDescendants[NodeIter->first].begin();

      oss << "[ ";
      for(; descentIter != this->NodeDescendants[NodeIter->first].end();
          ++descentIter)
        {
        std::string hcode = *descentIter;
        oss << this->GetNodeIndex(hcode) << " ";
        } // END for all descendants
      oss << "] ";
      if( this->NodeDescendants[NodeIter->first].size() > 1)
        {
        oss << "(**SPLIT**)";
        }
      if(NodeIter->second.HaloType == ZOMBIEHALO)
        {
        oss << "(**ZOMBIE**)";
        }
      }
    else
      {
      oss << "-1 (FINAL TIMESTEP)";
      }
    oss << std::endl;
    } // END for all nodes


  return( oss.str() );
}


} /* namespace cosmotk */
