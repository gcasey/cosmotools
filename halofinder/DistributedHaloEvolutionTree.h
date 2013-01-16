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
          int weight=1);

  /**
   * @brief Checks if the tree is empty.
   * @return true iff the tree is empty, else, false.
   */
  bool IsEmpty();

  /**
   * @brief Computes the number of nodes.
   * @return N the number of nodes.
   */
  int GetNumberOfNodes();

  /**
   * @brief Computes the number of edges.
   * @return N the number of edges.
   */
  int GetNunmberOfEdges();

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

  std::string FileName;

  MPI_Comm Communicator;

  /**
   * @brief Computes the range of IDs for this process via a prefix-sum
   * @param range the range (out)
   */
  void GetNodeRangeForProcess(ID_T range[2]);


private:
  DISABLE_COPY_AND_ASSIGNMENT(DistributedHaloEvolutionTree);

};

} /* namespace cosmotk */
#endif /* DISTRIBUTEDHALOEVOLUTIONTREE_H_ */
