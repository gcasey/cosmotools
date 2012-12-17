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
      Halo **halos, const int N)
{
  for( int halo=0; halo < N; ++halo )
    {
    this->Nodes[ halos[ halo ]->GetHashCode() ] = *(halos[ halo ]);
    } // END for all halos
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::CreateEdge(
      const std::string halo1,
      const std::string halo2,
      int w)
{
  assert("pre: node halo1 does not exist" && this->HasNode(halo1) );
  assert("pre: node halo2 does not exist" && this->HasNode(halo2) );
  this->Edges.push_back( halo1 );
  this->Edges.push_back( halo2 );
  this->EdgeWeights.push_back( w );
}

//------------------------------------------------------------------------------
bool DistributedHaloEvolutionTree::IsEmpty()
{
  return( this->Nodes.empty() );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfNodes()
{
  return( this->Nodes.size() );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNunmberOfEdges()
{
  return( this->Edges.size()/2 );
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
  assert("pre: NULL communicator!" && (this->Communicator != MPI_COMM_NULL));
  // TODO: need to implement this
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::WriteTree()
{
  assert("pre: NULL communicator!" && (this->Communicator != MPI_COMM_NULL));
  this->RelabelTreeNodes();
  // TODO: implement this
}



} /* namespace cosmotk */
