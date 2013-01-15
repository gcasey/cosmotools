#include "DistributedHaloEvolutionTree.h"

#include "Halo.h"

#include <cassert>

namespace cosmotk
{

DistributedHaloEvolutionTree::DistributedHaloEvolutionTree()
{
  this->Communicator = MPI_COMM_NULL;
  this->FileName = "HaloEvolutionTree.dat";
}

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::~DistributedHaloEvolutionTree()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::AppendNodes(
      Halo *halos, const int N)
{
  for( int halo=0; halo < N; ++halo )
    {
    this->Nodes[ halos[ halo ].GetHashCode() ] = halos[ halo ];

    // Remove the particle Ids of the object
    this->Nodes[ halos[ halo ].GetHashCode() ].ParticleIds.clear();
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
void DistributedHaloEvolutionTree::GetNodeRangeForProcess(int r[2])
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
  // TODO: Write tree in this->FileName
}



} /* namespace cosmotk */
