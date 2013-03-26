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
#include <set>    // For STL set

// MPI
#include <mpi.h>

// DIY
#include "diy.h" // For DIY_Datatype

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
   * @brief Returns the number of zombie nodes at the given timestep
   * @param tstep the timestep in query.
   * @return N the number of zombie nodes in the tree at the given timestep.
   */
  int GetNumberOfZombieNodes(const int tstep);

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
   * @brief Clears out all the data in this tree
   * @post this->IsEmpty()==true
   */
  void Clear();

  /**
   * @brief Barrier synchronization among all processes
   */
  void Barrier() { MPI_Barrier(this->Communicator); }

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

protected:
  std::map< std::string, Halo > Nodes; // List of nodes in the halo
  std::map< std::string, ID_T > Node2UniqueIdx; // Maps halo nodes to a idx
  std::vector< std::string >    Edges; // List of edges (strided by 2)
  std::vector< int > EdgeWeights;      // Weights associated with edges
  std::vector< int > EdgeEvents;       // edge events

  std::map<int,int> NodeCounter; // Counts number of nodes at a given timestep
  std::map<int,int> ZombieCounter; // Counts number of zombies at a given tstep
  std::string FileName;

  std::map< std::string, std::set<std::string> > NodeDescendants;

  MPI_Comm Communicator;

  int IOFormat;

  /**
   * @brief Updates the node counter at the given timestep
   * @param tstep the timestep to update.
   */
  void UpdateNodeCounter(const int tstep);

  /**
   * @brief Update the zombie counter at the given timestep.
   * @param tstep the timestep to update.
   */
  void UpdateZombieCounter(const int tstep);

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
