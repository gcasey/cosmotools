#include "DistributedHaloEvolutionTree.h"

#include "Halo.h"

#include <cassert>

namespace cosmotk
{

DistributedHaloEvolutionTree::DistributedHaloEvolutionTree()
{
  this->Communicator = MPI_COMM_NULL;
}

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::~DistributedHaloEvolutionTree()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::AppendNodes(
      std::vector< Halo >& halos)
{
  for( int halo=0; halo < static_cast<int>(halos.size()); ++halo )
    {
    this->Nodes[ halos[ halo ].GetHashCode() ] = halos[ halo ];
    } // END for all halos
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::CreateEdge(
      const std::string halo1,
      const std::string halo2  )
{
  assert("pre: node halo1 does not exist" && this->HasNode(halo1) );
  assert("pre: node halo2 does not exist" && this->HasNode(halo2) );
  this->Edges.push_back( halo1 );
  this->Edges.push_back( halo2 );
}

//------------------------------------------------------------------------------
bool DistributedHaloEvolutionTree::IsEmpty()
{
  return( this->Nodes.empty() );
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
void DistributedHaloEvolutionTree::RelabelTreeNodes()
{
  // TODO: need to implement this
}

} /* namespace cosmotk */
