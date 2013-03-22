#include "DistributedHaloEvolutionTree.h"

// CosmologyTools includes
#include "GenericIO.h"
#include "Halo.h"

// STL includes
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <set>

//------------------------------------------------------------------------------
// DIY callback routines
//------------------------------------------------------------------------------
void* create_write_type(void* block, int did, int lid, DIY_Datatype *dtype)
{
  DIYTree *treePtr = static_cast<DIYTree*>(block);
  cosmotk::DistributedHaloEvolutionTree::CreateDIYTreeType(treePtr,dtype);
  return DIY_BOTTOM;
}

void* create_read_type(
          int did, int lid, int *hdr, void **base_addr, DIY_Datatype *dtype)
{
  DIYTree *treePtr       = new DIYTree();
  treePtr->NumberOfNodes = hdr[0];
  treePtr->NumberOfEdges = hdr[1];
  treePtr->TreeNodes = new DIYTreeNodeType[treePtr->NumberOfNodes];
  treePtr->TreeEdges = new DIYTreeEdgeType[treePtr->NumberOfEdges];
  cosmotk::DistributedHaloEvolutionTree::CreateDIYTreeType(treePtr,dtype);
  *base_addr = DIY_BOTTOM;
  return treePtr;
}

//------------------------------------------------------------------------------

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
  this->IOFormat = MergerTreeFileFormat::GENERIC_IO;
}

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::~DistributedHaloEvolutionTree()
{
  this->Clear();
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::GetDIYTree(DIYTree* tree)
{
  assert("pre: DIYTree pointer is NULL!" && (tree != NULL));

  // STEP 0: Allocate data-structures
  tree->NumberOfNodes = this->GetNumberOfNodes();
  tree->NumberOfEdges = this->GetNumberOfEdges();
  tree->TreeNodes     = new DIYTreeNodeType[tree->NumberOfNodes];
  tree->TreeEdges     = new DIYTreeEdgeType[tree->NumberOfEdges];

  // STEP 1: Fill nodes
  std::map<std::string,Halo>::iterator iter = this->Nodes.begin();
  for(int idx=0; iter != this->Nodes.end(); ++iter, ++idx)
    {
    std::string hcode                 = iter->first;
    tree->TreeNodes[idx].HaloTag      = iter->second.Tag;
    tree->TreeNodes[idx].TimeStep     = iter->second.TimeStep;
    tree->TreeNodes[idx].RedShift     = iter->second.Redshift;
    tree->TreeNodes[idx].UniqueNodeID = this->GetNodeIndex( hcode );
    for( int i=0; i < 3; ++i )
      {
      tree->TreeNodes[idx].HaloCenter[i]   = iter->second.Center[i];
      tree->TreeNodes[idx].HaloVelocity[i] = iter->second.AverageVelocity[i];
      }

    } // END for all nodes

  // STEP 2: Fill edges
  for(int edgeIdx=0; edgeIdx < this->GetNumberOfEdges(); ++edgeIdx)
    {
    ID_T node1 = this->GetNodeIndex( this->Edges[edgeIdx*2]   );
    ID_T node2 = this->GetNodeIndex( this->Edges[edgeIdx*2+1] );
    tree->TreeEdges[edgeIdx].EndNodes[0] = node1;
    tree->TreeEdges[edgeIdx].EndNodes[1] = node2;
    tree->TreeEdges[edgeIdx].EdgeWeight  = this->EdgeWeights[edgeIdx];
    tree->TreeEdges[edgeIdx].EdgeEvent   = this->EdgeEvents[edgeIdx];
    } // END for all edges
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

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::WriteTree()
{
  assert("pre: NULL communicator!" && (this->Communicator != MPI_COMM_NULL));
  this->RelabelTreeNodes();

  switch(this->IOFormat)
    {
    case MergerTreeFileFormat::DIY:
      this->WriteWithDIY();
      break;
    case MergerTreeFileFormat::GENERIC_IO:
      this->WriteWithGenericIO();
      break;
    default:
      std::cerr << "ERROR: Undefined IO format!" << std::endl;
      assert(false);
    }
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::ReadTree()
{
  switch(this->IOFormat)
    {
    case MergerTreeFileFormat::DIY:
      this->ReadWithDIY();
      break;
    case MergerTreeFileFormat::GENERIC_IO:
      this->ReadWithGenericIO();
      break;
    default:
      std::cerr << "ERROR: Undefined IO format!" << std::endl;
      assert(false);
    }
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::ReadWithDIY()
{
  // STEP 0: Open file to read
  int swap     = 0;
  int compress = 0;
  int numblocks =
      DIY_Read_open_all(
          0,const_cast<char*>(this->FileName.c_str()),swap,compress);

  // STEP 1: Allocate headers for blocks
  void **pmblocks;
  int **hdrs = new int*[numblocks]; // How do we know the number of blocks a priori?
  for(int i=0; i < numblocks; ++i)
    {
    hdrs[i] = new int[DIY_MAX_HDR_ELEMENTS];
    }

  // STEP 2: Read blocks
  DIY_Read_blocks_all(0,&pmblocks,hdrs,&create_read_type);

  // STEP 3: Unpack
  // TODO: implement this

  // STEP 4: Close file
  DIY_Read_close_all(0);
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::WriteWithDIY()
{
  // STEP 0: Get DIYTree representation
  DIYTree *mtree = new DIYTree;
  this->GetDIYTree( mtree );

  // STEP 1: Pack DIY tree
  int nblocks = 1;
  int **hdrs  = new int*[nblocks];
  void **pmblocks;
  pmblocks = new void*[nblocks];
  for( int blkIdx=0; blkIdx < nblocks; ++blkIdx)
    {
    hdrs[blkIdx]     = new int[DIY_MAX_HDR_ELEMENTS];
    hdrs[blkIdx][0]  = mtree->NumberOfNodes;
    hdrs[blkIdx][1]  = mtree->NumberOfEdges;
    pmblocks[blkIdx] = mtree;
    } // END for all blocks

  // STEP 2: DIY write
  int compress = 0;
  DIY_Write_open_all(0,const_cast<char*>(this->FileName.c_str()),compress);
  DIY_Write_blocks_all(0,pmblocks,nblocks,hdrs,&create_write_type);
  DIY_Write_close_all(0);

  // STEP 3: Clean up
  delete mtree;
  for(int blkIdx=0; blkIdx < nblocks; ++blkIdx)
    {
    delete [] hdrs[blkIdx];
    } // END for all blocks
  delete [] hdrs;
  delete [] pmblocks;
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::GetFlatTreeNodeArrays(
      ID_T *nodeIdx, ID_T *haloTags, ID_T *tstep, REAL *redShift,
      POSVEL_T *center_x, POSVEL_T *center_y, POSVEL_T *center_z,
      POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz)
{
  assert("pre: nodeIdx array is NULL!"  && (nodeIdx != NULL) );
  assert("pre: haloTags array is NULL!" && (haloTags != NULL) );
  assert("pre: redShift array is NULL!" && (redShift != NULL) );
  assert("pre: center_x array is NULL!" && (center_x != NULL) );
  assert("pre: center_y array is NULL!" && (center_y != NULL) );
  assert("pre: center_z array is NULL!" && (center_z != NULL) );
  assert("pre: vx array is NULL!" && (vx != NULL) );
  assert("pre: vy array is NULL!" && (vy != NULL) );
  assert("pre: vz array is NULL!" && (vz != NULL) );

  int idx = 0;
  std::map<std::string,Halo>::iterator NodesIter = this->Nodes.begin();
  for(; NodesIter != this->Nodes.end(); ++NodesIter, ++idx)
    {
    std::string hashCode = NodesIter->first;
    nodeIdx[ idx ]  = this->GetNodeIndex( hashCode );
    haloTags[ idx ] = NodesIter->second.Tag;
    tstep[ idx ]    = NodesIter->second.TimeStep;
    redShift[ idx ] = NodesIter->second.Redshift;
    center_x[ idx ] = NodesIter->second.Center[0];
    center_y[ idx ] = NodesIter->second.Center[1];
    center_z[ idx ] = NodesIter->second.Center[2];
    vx[ idx ] = NodesIter->second.AverageVelocity[0];
    vy[ idx ] = NodesIter->second.AverageVelocity[1];
    vz[ idx ] = NodesIter->second.AverageVelocity[2];
    } // END for all nodes
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::GetFlatTreeEdgeArrays(
        ID_T *startNodeIdx, ID_T *endNodeIdx )
{
  assert("pre: startNodeIdx array is NULL!" && (startNodeIdx != NULL) );
  assert("pre: endNodeIdx array is NULL!" && (endNodeIdx != NULL) );

  for( int edgeIdx=0; edgeIdx < this->GetNumberOfEdges(); ++edgeIdx )
    {
    startNodeIdx[edgeIdx] = this->GetNodeIndex(this->Edges[edgeIdx*2]);
    endNodeIdx[edgeIdx]   = this->GetNodeIndex(this->Edges[edgeIdx*2+1]);
    } // END for all edges
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::ReadWithGenericIO()
{
//  // STEP 0: Clear data
//  this->Clear();
//
//  // STEP 1: Initialize Nodes and Edge reader
//  std::ostringstream oss;
//  oss << this->FileName << ".nodes";
//  std::string nodesFile = oss.str();
//  cosmotk::GenericIO GIONodeReader(this->Communicator,nodesFile);
//
//  oss.clear(); oss.str("");
//  oss << this->FileName << ".edges";
//  std::string edgesFile = oss.str();
//  cosmotk::GenericIO GIOEdgesReader(this->Communicator,edgesFile);
//
//  // STEP 2: Read headers
//  GIONodeReader.openAndReadHeader(false);
//  GIOEdgesReader.openAndReadHeader(false);
//  assert( "pre: nodes file and edges file rank mismatch!" &&
//   GIONodeReader.readNRanks() == GIOEdgesReader.readNRanks() );
//
//  // STEP 3: Assign blocks to processes
//  int numBlocks = GIONodeReader.readNRanks();
//  std::vector<int> assignedBlocks;
//  this->RoundRobinAssignment(numBlocks,assignedBlocks);
//
//  // STEP 4: Read in the data for rank in file
//  for(unsigned int i=0; i < assignedBlocks.size(); ++i)
//    {
//    int block = assignedBlocks[i];
//    this->ReadBlock(block,&GIONodeReader,&GIOEdgesReader);
//    } // END for all ranks
//
//  GIONodeReader.close();
//  GIOEdgesReader.close();
}

//------------------------------------------------------------------------------
//void DistributedHaloEvolutionTree::ReadBlock(
//        int block,
//        cosmotk::GenericIO *nodesReader,
//        cosmotk::GenericIO *edgesReader)
//{
//  assert("pre: Nodes reader is NULL!" && (nodesReader != NULL) );
//  assert("pre: Nodes reader is NULL!" && (edgesReader != NULL) );
//
//  // STEP 0: Get number of nodes
//  int numNodes = nodesReader->readNumElems(block);
//
//  // STEP 1: Allocate temporary arrays
//  // TODO: Should I allocate extra space here? Not sure how?
//  ID_T *nodeIdx      = new ID_T[numNodes];
//  ID_T *haloTags     = new ID_T[numNodes];
//  ID_T *tstep        = new ID_T[numNodes];
//  REAL *redShift     = new REAL[numNodes];
//  POSVEL_T *center_x = new POSVEL_T[numNodes];
//  POSVEL_T *center_y = new POSVEL_T[numNodes];
//  POSVEL_T *center_z = new POSVEL_T[numNodes];
//  POSVEL_T *vx = new POSVEL_T[numNodes];
//  POSVEL_T *vy = new POSVEL_T[numNodes];
//  POSVEL_T *vz = new POSVEL_T[numNodes];
//
//  // STEP 2: Add variables to be read
//  nodesReader->addVariable("UniqueIndex", nodeIdx, true);
//  nodesReader->addVariable("HaloTags", haloTags, true);
//  nodesReader->addVariable("TimeStep", tstep, true);
//  nodesReader->addVariable("Red-Shifts", redShift, true);
//  nodesReader->addVariable("Center-X", center_x, true);
//  nodesReader->addVariable("Center-Y", center_y, true);
//  nodesReader->addVariable("Center-Z", center_z, true);
//  nodesReader->addVariable("Velocity-X", vx, true);
//  nodesReader->addVariable("Velocity-Y", vy, true);
//  nodesReader->addVariable("Velocity-Z", vz, true);
//
//  // STEP 3: Do the I/O
//  nodesReader->readData(block,false,false);
//
//  // STEP 4: Unpack the nodes
//  std::map< ID_T, std::string > index2hashcode;
//  cosmotk::Halo h;
//  for( int idx=0; idx < numNodes; ++idx)
//    {
//    h.Tag = haloTags[ idx ];
//    h.TimeStep = tstep[ idx ];
//    h.Redshift = redShift[ idx ];
//    h.Center[0] = center_x[ idx ];
//    h.Center[1] = center_y[ idx ];
//    h.Center[2] = center_z[ idx ];
//    h.AverageVelocity[0] = vx[ idx ];
//    h.AverageVelocity[1] = vy[ idx ];
//    h.AverageVelocity[2] = vz[ idx ];
//    index2hashcode[ nodeIdx[idx] ] = h.GetHashCode();
//    this->Node2UniqueIdx[ h.GetHashCode() ] = nodeIdx[idx];
//    this->Nodes[ h.GetHashCode() ] = h;
//    } // END for all nodes
//
//  // STEP 5: Clean Up temporary flat arrays
//  delete [] nodeIdx;
//  delete [] haloTags;
//  delete [] tstep;
//  delete [] redShift;
//  delete [] center_x;
//  delete [] center_y;
//  delete [] center_z;
//  delete [] vx;
//  delete [] vy;
//  delete [] vz;
//
//  // STEP 6: Get number of edges
//  int numEdges = edgesReader->readNumElems(block);
//
//  // STEP 7: Allocate and acquire flat edge arrays
//  // TODO: Should I allocate extra space here? Not sure how?
//  ID_T *startNodeIdx = new ID_T[numEdges];
//  ID_T *endNodeIdx = new ID_T[numEdges];
//  int *edgeWeights = new int[numEdges];
//  int *edgeEvents  = new int[numEdges];
//
//
//  // STEP 8: Add variables to be read
//  edgesReader->addVariable("StartNode", startNodeIdx,true);
//  edgesReader->addVariable("EndNode", endNodeIdx,true);
//  edgesReader->addVariable("Weights", edgeWeights,true);
//  edgesReader->addVariable("Events", edgeEvents,true);
//
//  // STEP 9: Do the I/O
//  edgesReader->readData(block,false,false);
//
//  // STEP 10: UnPack
//  for(int idx=0; idx < numEdges; ++idx)
//    {
//    assert( "ERROR: cannot find hash-code for edge start index" &&
//      index2hashcode.find(startNodeIdx[idx]) != index2hashcode.end());
//    assert( "ERROR: cannot find hash-code for edge end index" &&
//      index2hashcode.find(endNodeIdx[idx]) != index2hashcode.end());
//    this->Edges.push_back(index2hashcode[startNodeIdx[idx]]);
//    this->Edges.push_back(index2hashcode[endNodeIdx[idx]]);
//    this->EdgeWeights.push_back( edgeWeights[idx] );
//    this->EdgeEvents.push_back( edgeEvents[idx] );
//    } // END for all edges
//
//  // STEP 11: Clean up temporary flat arrays
//  delete [] startNodeIdx;
//  delete [] endNodeIdx;
//  delete [] edgeWeights;
//  delete [] edgeEvents;
//
//  // STEP 12: Clear variables since the readers are called iteratively
//  nodesReader->clearVariables();
//  edgesReader->clearVariables();
//}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::WriteWithGenericIO()
{
  // STEP 0: Initialize Nodes writer
  std::ostringstream oss;
  oss << this->FileName << ".nodes";
  cosmotk::GenericIO GIONodeWriter(this->Communicator,oss.str());
  GIONodeWriter.setNumElems(this->GetNumberOfNodes());

  // STEP 1: Allocate and acquire flat node arrays
  ID_T *nodeIdx      = new ID_T[this->GetNumberOfNodes()];
  ID_T *haloTags     = new ID_T[this->GetNumberOfNodes()];
  ID_T *tstep        = new ID_T[this->GetNumberOfNodes()];
  REAL *redShift     = new REAL[this->GetNumberOfNodes()];
  POSVEL_T *center_x = new POSVEL_T[this->GetNumberOfNodes()];
  POSVEL_T *center_y = new POSVEL_T[this->GetNumberOfNodes()];
  POSVEL_T *center_z = new POSVEL_T[this->GetNumberOfNodes()];
  POSVEL_T *vx = new POSVEL_T[this->GetNumberOfNodes()];
  POSVEL_T *vy = new POSVEL_T[this->GetNumberOfNodes()];
  POSVEL_T *vz = new POSVEL_T[this->GetNumberOfNodes()];
  this->GetFlatTreeNodeArrays(
      nodeIdx,haloTags,tstep, redShift,
      center_x,center_y,center_z,
      vx,vy,vz);

  // STEP 2: Add variables to be written
  GIONodeWriter.addVariable("UniqueIndex", nodeIdx);
  GIONodeWriter.addVariable("HaloTags", haloTags);
  GIONodeWriter.addVariable("TimeStep", tstep);
  GIONodeWriter.addVariable("Red-Shifts", redShift);
  GIONodeWriter.addVariable("Center-X", center_x);
  GIONodeWriter.addVariable("Center-Y", center_y);
  GIONodeWriter.addVariable("Center-Z", center_z);
  GIONodeWriter.addVariable("Velocity-X", vx);
  GIONodeWriter.addVariable("Velocity-Y", vy);
  GIONodeWriter.addVariable("Velocity-Z", vz);

  // STEP 3: Do the I/O
  GIONodeWriter.write();

  // STEP 4: Clean Up temporary flat arrays
  delete [] nodeIdx;
  delete [] haloTags;
  delete [] tstep;
  delete [] redShift;
  delete [] center_x;
  delete [] center_y;
  delete [] center_z;
  delete [] vx;
  delete [] vy;
  delete [] vz;

  // STEP 5: Initialize edges writer
  oss.clear(); oss.str("");
  oss << this->FileName << ".edges";
  cosmotk::GenericIO GIOEdgeWriter(this->Communicator,oss.str());
  GIOEdgeWriter.setNumElems(this->GetNumberOfEdges());

  // STEP 6: Allocate and acquire flat edge arrays
  ID_T *startNodeIdx = new ID_T[this->GetNumberOfEdges()];
  ID_T *endNodeIdx = new ID_T[this->GetNumberOfEdges()];
  this->GetFlatTreeEdgeArrays(startNodeIdx, endNodeIdx);

  // STEP 7: Add variables to be written
  GIOEdgeWriter.addVariable("StartNode", startNodeIdx);
  GIOEdgeWriter.addVariable("EndNode", endNodeIdx);
  GIOEdgeWriter.addVariable("Weights",&this->EdgeWeights[0]);
  GIOEdgeWriter.addVariable("Events", &this->EdgeEvents[0]);

  // STEP 8: Do the I/O
  GIOEdgeWriter.write();

  // STEP 9: Clean up temporary flat arrays
  delete [] startNodeIdx;
  delete [] endNodeIdx;
}

} /* namespace cosmotk */
