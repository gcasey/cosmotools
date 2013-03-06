/**
 * @brief A distributed halo evolution tree data-structure used to capture
 * the evolution of halos during the course of a simulation.
 *
 * A DistributedHaloEvolutionTree data-structure is essentially a distributed
 * directed graph. Each node in the tree corresponds to a halo and the edges
 * between nodes indicate parent/child relationship. The data-structure is
 * designed s.t. it can be constructed incrementally.
 */
#ifndef DISTRIBUTEDHALOEVOLUTIONTREE_H_
#define DISTRIBUTEDHALOEVOLUTIONTREE_H_

#include "CosmologyToolsMacros.h"
#include "MergerTreeEvent.h"
#include "MergerTreeFileFormat.h"

// C++ STL includes
#include <vector> // For STL vector
#include <map>    // For STL map

// MPI
#include <mpi.h>

// DIY
#include "diy.h" // For DIY_Datatype

/**
 * @struct DIYTreeNodeType
 * @brief Used to store tree nodes.
 */
struct DIYTreeNodeType {
  ID_T UniqueNodeID;
  ID_T HaloTag;
  ID_T TimeStep;
  REAL RedShift;
  POSVEL_T HaloCenter[3];
  POSVEL_T HaloVelocity[3];
};

/**
 * @struct DIYTreeEdgeType
 * @brief Used to store tree edges.
 */
struct DIYTreeEdgeType {
  ID_T EndNodes[2];
  int EdgeWeight;
  int EdgeEvent;
};

/**
 * @struct DIYTree
 * @brief Stores pointers to tree nodes and tree edges.
 */
struct DIYTree {
  DIYTreeNodeType* TreeNodes;

  // NOTE: Metadata, these fields are skipped when creating a DIY data type,
  // hence, they are injected in the middle of the struct, s.t. MPI will compute
  // the extents of the datatype correctly.
  int NumberOfNodes;
  int NumberOfEdges;

  DIYTreeEdgeType* TreeEdges;
};

namespace cosmotk
{

// Forward declarations
class Halo;
class GenericIO;

class DistributedHaloEvolutionTree
{
public:
  DistributedHaloEvolutionTree();
  virtual ~DistributedHaloEvolutionTree();

  /**
   * @brief Get/Set macro for MPI communicator.
   */
  GetNSetMacro(Communicator,MPI_Comm);

  /**
   * @brief Get/Set macro for FileName
   * @note This is the filename where the merger will be written out.
   */
  GetNSetMacro(FileName,std::string);

  /**
   * @brief Get/Set macro for the MergerTree format
   * @note Default is set to GENERIC_IO
   * @see enum MergerTreeFormat
   */
  GetNSetMacro(IOFormat,int);

  /**
   * @brief Inserts the given halo as a node in the tree.
   * @param halo the halo to add
   */
  void InsertNode(Halo &halo);

  /**
   * @brief Adds the given halos as nodes to this instance of the tree
   * @param halos the set of halos to be added in the tree.
   * @note Each tree node is identified by the halo hash-code.
   * @see Halo::GetHashCode()
   */
  void AppendNodes( Halo* halos, const int N );

  /**
   * @brief Creates an edge from the tree node corresponding to hashCode1
   * to the tree node corresponding to hashCode2.
   * @param hashCode1 the hashCode of the 1st tree node.
   * @param hashCode2 the hashCode of the 2nd tree node.
   * @pre this->HasNode( hashCode1 ) && this->HasNode( hashCode2 ) is true.
   */
  void CreateEdge(
          std::string hashCode1,
          std::string hashCode2,
          int weight=1,
          int event=MergerTreeEvent::UNDEFINED);

  /**
   * @brief Checks if the tree is empty.
   * @return true iff the tree is empty, else, false.
   */
  bool IsEmpty();

  /**
   * @brief Returns the number of nodes at the given timestep.
   * @param tstep the timestep in query.
   * @return N the number of nodes in the tree at the given timestep.
   */
  int GetNumberOfNodes(const int tstep);

  /**
   * @brief Computes the number of nodes.
   * @return N the number of nodes.
   */
  int GetNumberOfNodes();

  /**
   * @brief Computes the number of edges.
   * @return N the number of edges.
   */
  int GetNumberOfEdges();

  /**
   * @brief Checks if the halo corresponding to the given hashCode exists
   * (in this process) on the tree.
   * @param hashCode the hash code corresponding to the halo in query.
   * @return true iff the halo exists in the tree, else false.
   */
  bool HasNode( std::string hashCode1 );

  /**
   * @brief Relabels the tree and dumps it into a file.
   */
  void WriteTree();

  /**
   * @brief Reads the tree from a file.
   */
  void ReadTree();

  /**
   * @brief Clears out all the data in this tree
   * @post this->IsEmpty()==true
   */
  void Clear();

  /**
   * @brief Barrier synchronization among all processes
   */
  void Barrier() { MPI_Barrier(this->Communicator); }

  /**
   * @brief Get DIYTree instance for the halo-evolution tree
   * @param tree pointer to the DIYTree struct (in/out)
   */
  void GetDIYTree(DIYTree* tree);

  /**
   * @brief Returns the unique ID for a given node
   * @param hashCode the hashcode for the node
   * @return idx the ID of the node
   */
  ID_T GetNodeIndex(std::string hashCode);

  /**
   * @brief Relabels s.t. all nodes are uniquely identified by an integer ID
   * across all processes.
   */
  void RelabelTreeNodes();

  /**
   * @brief Creates an ASCII string representation of the given tree
   * @return a string representing the data in the tree
   * @note Used mainly for debugging.
   */
  std::string ToString();

  /**
   * @brief Registers a DIY data-type to represent a DIY tree.
   * @param myTree pointer to a DIY tree instance (in)
   * @param dtype pointer to the DIY data type (out)
   * @note This DIY datatype should only be used in the context of
   * DIY I/O since each type has variable sizes
   */
  static void CreateDIYTreeType(DIYTree *myTree, DIY_Datatype *dtype);

protected:
  std::map< std::string, Halo > Nodes; // List of nodes in the halo
  std::map< std::string, ID_T > Node2UniqueIdx; // Maps halo nodes to a idx
  std::vector< std::string >    Edges; // List of edges (strided by 2)
  std::vector< int > EdgeWeights;      // Weights associated with edges
  std::vector< int > EdgeEvents;       // edge events

  std::map<int,int> NodeCounter; // Counts number of nodes at a given timestep
  std::string FileName;

  MPI_Comm Communicator;

  int IOFormat;

  /**
   * @brief Updates the node counter at the given timestep
   * @param tstep the timestep to update.
   */
  void UpdateNodeCounter(const int tstep);

  /**
   * @brief Gets the variables associated with each tree edge in flat arrays.
   * @param startNodeIdx list of start node indices
   * @param endNodeIdx list of end node indices
   * @param weights list of edge weights
   * @param eventType list of edge event types
   */
  void GetFlatTreeEdgeArrays(ID_T *startNodeIdx, ID_T *endNodeIdx);

  /**
   * @brief Gets the variables associated with each tree node in flat arrays.
   * @param nodeIds the unique index (computed) of node (i.e.,halo) in the tree
   * @param haloTags the original halo tag (computed from the halo-finder)
   * @param redShift the red-shift of the associate halo
   * @param center_x the x-coordinate of the halo center
   * @param center_y the y-coordinate of the halo center
   * @param center_z the z-coordinate of the halo center
   * @param vx the x-coordinate of the halo velocity vector
   * @param vy the y-coordinate of the halo velocity vector
   * @param vz the z-coordinate of the halo velocity vector
   * @pre All arrays that are passed in must be pre-allocated
   * @note This method is used to create the flat arrays needed for the
   * GenericIO
   */
  void GetFlatTreeNodeArrays(
      ID_T *nodeIds,
      ID_T *haloTags,
      ID_T *timestep,
      REAL *redShift,
      POSVEL_T *center_x, POSVEL_T *center_y, POSVEL_T *center_z,
      POSVEL_T *vx, POSVEL_T *vy, POSVEL_T *vz);

  /**
   * @brief Computes the range of IDs for this process via a prefix-sum
   * @param range the range (out)
   */
  void GetNodeRangeForProcess(ID_T range[2]);

  /**
   * @brief Reads data in DIY format
   */
  void ReadWithDIY();

  /**
   * @brief Writes data in DIY format
   */
  void WriteWithDIY();

  /**
   * @brief Reads data in GenericIO format
   */
  void ReadWithGenericIO();

  /**
   * @brief Reads in the given block
   * @param block the block to read
   * @param nodesReader pointer to the nodes reader
   * @param edgesReader pointer to the edges reader
   * @pre nodesReader != NULL
   * @pre edgesReader != NULL
   * @note Called from ReadWithGenericIO to read in a
   * particular block to handle the case where each process
   * reads more than one block.
   */
//  void ReadBlock(
//      int block, cosmotk::GenericIO *nodesReader,
//      cosmotk::GenericIO *edgesReader);

  /**
   * @brief Writes data in GenericIO format
   */
  void WriteWithGenericIO();

private:
  DISABLE_COPY_AND_ASSIGNMENT(DistributedHaloEvolutionTree);

};

} /* namespace cosmotk */
#endif /* DISTRIBUTEDHALOEVOLUTIONTREE_H_ */
