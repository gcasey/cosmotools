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
}

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::~DistributedHaloEvolutionTree()
{
  this->Clear();
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::InsertNode(
      const HaloInfo &halo, const unsigned char mask)
{
//  this->Nodes[ halo.GetHashCode() ] = halo;
//  this->Nodes[ halo.GetHashCode() ].ParticleIds.clear();
//  if( halo.HaloType != GHOSTHALO )
//    {
//    this->UpdateNodeCounter( halo.TimeStep );
//    }
//
//  if( halo.HaloType == ZOMBIEHALO )
//    {
//    this->UpdateZombieCounter( halo.TimeStep );
//    }
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::LinkHalos(
      const std::string halo1,
      const std::string halo2)
{
//  assert("pre: node halo1 does not exist" && this->HasNode(halo1) );
//  assert("pre: node halo2 does not exist" && this->HasNode(halo2) );
//  this->Edges.push_back( halo1 );
//  this->Edges.push_back( halo2 );
//  this->EdgeWeights.push_back( w );
//  this->EdgeEvents.push_back( e );
//
//  if( this->NodeDescendants.find(halo1) != this->NodeDescendants.end() )
//    {
//    this->NodeDescendants[ halo1 ].insert( halo2 );
//    }
//  else
//    {
//    std::set< std::string > descendants;
//    descendants.insert( halo2 );
//    this->NodeDescendants[ halo1 ] = descendants;
//    }
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::Clear()
{
//  this->Nodes.clear();
//  this->Node2UniqueIdx.clear();
//  this->Edges.clear();
//  this->EdgeWeights.clear();
//  this->EdgeEvents.clear();
}

//------------------------------------------------------------------------------
bool DistributedHaloEvolutionTree::IsEmpty()
{
  // TODO: implement this
  return false;
//  return( this->Nodes.empty() );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfNodes(const int tstep)
{
  return( -1 );
//  int N = 0;
//  if(this->NodeCounter.find(tstep) != this->NodeCounter.end() )
//    {
//    N = this->NodeCounter[ tstep ];
//    }
//  else
//    {
//    std::cerr << "WARNING: No tree nodes for requrested time-step!\n";
//    }
//  return( N );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfZombieNodes(
    const int tstep)
{
  return( -1 );
//  int N = 0;
//  if( this->ZombieCounter.find(tstep) != this->ZombieCounter.end() )
//    {
//    N = this->ZombieCounter[ tstep ];
//    }
//  return( N );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfNodes()
{
  return( -1 );
//  return( this->Nodes.size() );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfEdges()
{
  return( -1 );
//  return( this->EdgeWeights.size() );
}

//------------------------------------------------------------------------------
bool DistributedHaloEvolutionTree::HasNode(std::string hashCode)
{
  return false;
//  if( this->Nodes.find(hashCode) == this->Nodes.end() )
//    {
//    return false;
//    }
//  return true;
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
ID_T DistributedHaloEvolutionTree::GetNodeIndex(std::string hashCode)
{
  return( -1 );
//  assert( "pre: Cannot find Unique index for node!" &&
//      this->Node2UniqueIdx.find(hashCode)!=this->Node2UniqueIdx.end());
//  return(this->Node2UniqueIdx[hashCode]);
}

//------------------------------------------------------------------------------
std::string DistributedHaloEvolutionTree::ToString()
{
//  this->RelabelTreeNodes();
//
//  std::ostringstream oss;
//  oss << "NUMNODES " << this->GetNumberOfNodes() << std::endl;
//  std::map<std::string, Halo>::iterator NodeIter = this->Nodes.begin();
//  oss << "ID\t";
//  oss << "TIMESTEP\t";
//  oss << "TAG\t";
//  oss << "REDSHIFT\t";
//  oss << "CENTER_X\t";
//  oss << "CENTER_Y\t";
//  oss << "CENTER_Z\t";
//  oss << "MEAN_CENTER_X\t";
//  oss << "MEAN_CENTER_Y\t";
//  oss << "MEAN_CENTER_Z\t";
//  oss << "V_X\t";
//  oss << "V_Y\t";
//  oss << "V_Z\t";
//  oss << "MASS\t";
//  oss << "DESCENDANTS";
//  oss << std::endl;
//  for(;NodeIter != this->Nodes.end(); ++NodeIter)
//    {
//    oss << std::scientific
//        << std::setprecision(std::numeric_limits<POSVEL_T>::digits10);
//    oss << this->GetNodeIndex(NodeIter->first) << "\t";
//    oss << NodeIter->second.TimeStep  << "\t";
//    oss << NodeIter->second.Tag       << "\t";
//    oss << NodeIter->second.Redshift  << "\t";
//    oss << NodeIter->second.Center[0] << "\t";
//    oss << NodeIter->second.Center[1] << "\t";
//    oss << NodeIter->second.Center[2] << "\t";
//    oss << NodeIter->second.MeanCenter[0] << "\t";
//    oss << NodeIter->second.MeanCenter[1] << "\t";
//    oss << NodeIter->second.MeanCenter[2] << "\t";
//    oss << NodeIter->second.AverageVelocity[0] << "\t";
//    oss << NodeIter->second.AverageVelocity[1] << "\t";
//    oss << NodeIter->second.AverageVelocity[2] << "\t";
//    oss << NodeIter->second.HaloMass << "\t";
//
//    if( this->NodeDescendants.find(NodeIter->first) !=
//        this->NodeDescendants.end())
//      {
//      std::set< std::string >::iterator descentIter =
//          this->NodeDescendants[NodeIter->first].begin();
//
//      oss << "[ ";
//      for(; descentIter != this->NodeDescendants[NodeIter->first].end();
//          ++descentIter)
//        {
//        std::string hcode = *descentIter;
//        oss << this->GetNodeIndex(hcode) << " ";
//        } // END for all descendants
//      oss << "] ";
//      if( this->NodeDescendants[NodeIter->first].size() > 1)
//        {
//        oss << "(**SPLIT**)";
//        }
//      if(NodeIter->second.HaloType == ZOMBIEHALO)
//        {
//        oss << "(**ZOMBIE**)";
//        }
//      }
//    else
//      {
//      oss << "-1 (FINAL TIMESTEP)";
//      }
//    oss << std::endl;
//    } // END for all nodes
//
//
//  return( oss.str() );
}


} /* namespace cosmotk */
