#include "DistributedHaloEvolutionTree.h"

#include "Halo.h"

#include <cstdio>
#include <cstddef>
#include <cassert>

namespace cosmotk
{

void DistributedHaloEvolutionTree::CreateDIYTreeType(
    DIYTree *myTree,
    DIY_Datatype *dtype)
{
  assert("pre: DIYTree pointer is NULL!" && (myTree != NULL) );

  // STEP 0: Create the node type
  DIY_Datatype nodeType;
  struct map_block_t tree_node_map[6] = {
    {DIY_ID_T,     OFST, 1, offsetof(struct DIYTreeNodeType,UniqueNodeID)},
    {DIY_ID_T,     OFST, 1, offsetof(struct DIYTreeNodeType,HaloTag)},
    {DIY_ID_T,     OFST, 1, offsetof(struct DIYTreeNodeType,TimeStep)},
    {DIY_REAL_T,   OFST, 1, offsetof(struct DIYTreeNodeType,RedShift)},
    {DIY_POSVEL_T, OFST, 3, offsetof(struct DIYTreeNodeType,HaloCenter)},
    {DIY_POSVEL_T, OFST, 3, offsetof(struct DIYTreeNodeType,HaloVelocity)},
  };
 DIY_Create_struct_datatype(0, 6, tree_node_map, &nodeType);

 // STEP 1: Create the edge type
 DIY_Datatype edgeType;
 struct map_block_t edge_map[3] = {
   {DIY_ID_T, OFST, 2, offsetof(struct DIYTreeEdgeType,EndNodes)},
   {DIY_INT, OFST, 1, offsetof(struct DIYTreeEdgeType,EdgeWeight)},
   {DIY_INT, OFST, 1, offsetof(struct DIYTreeEdgeType,EdgeEvent)},
 };
 DIY_Create_struct_datatype(0, 3, edge_map, &edgeType);

 // STEP 2: Create the tree type
 struct map_block_t tree_map[2] = {
   {nodeType, ADDR, myTree->NumberOfNodes, DIY_Addr(myTree->TreeNodes)},
   {edgeType, ADDR, myTree->NumberOfEdges, DIY_Addr(myTree->TreeEdges)},
 };
 DIY_Create_struct_datatype(DIY_Addr(myTree), 2, tree_map, dtype);

 // STEP 3: Clean up
 DIY_Destroy_datatype(&nodeType);
 DIY_Destroy_datatype(&edgeType);
}

//------------------------------------------------------------------------------
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
void DistributedHaloEvolutionTree::GetDIYTree(DIYTree* tree)
{
  assert("pre: DIYTree pointer is NULL!" && (tree != NULL));

  // STEP 0: Allocate data-structures
  tree->NumberOfNodes = this->GetNumberOfNodes();
  tree->NumberOfEdges = this->GetNunmberOfEdges();
  tree->TreeNodes     = new DIYTreeNodeType[tree->NumberOfNodes];
  tree->TreeEdges     = new DIYTreeEdgeType[tree->NumberOfEdges];

  // STEP 1: Fill nodes
  std::map<std::string,Halo>::iterator iter = this->Nodes.begin();
  for(int idx=0; iter != this->Nodes.end(); ++iter, ++idx)
    {
    std::string hcode = iter->first;
    } // END for all nodes

  // STEP 2: Fill edges
  for(int edgeIdx=0; edgeIdx < this->GetNunmberOfEdges(); ++edgeIdx)
    {
    ID_T node1 = this->GetNodeIndex( this->Edges[edgeIdx*2]   );
    ID_T node2 = this->GetNodeIndex( this->Edges[edgeIdx*2+1] );
    tree->TreeEdges[edgeIdx].EndNodes[0] = node1;
    tree->TreeEdges[edgeIdx].EndNodes[1] = node2;
    tree->TreeEdges[edgeIdx].EdgeWeight  = this->EdgeWeights[edgeIdx];

    // TODO: Handle edge events correctly
    tree->TreeEdges[edgeIdx].EdgeEvent  = 0;
    } // END for all edges
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
void DistributedHaloEvolutionTree::WriteTree()
{
  assert("pre: NULL communicator!" && (this->Communicator != MPI_COMM_NULL));
  this->RelabelTreeNodes();
  // TODO: Write tree in this->FileName
}



} /* namespace cosmotk */
