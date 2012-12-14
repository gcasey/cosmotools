/**
 * @brief
 */
#ifndef DISTRIBUTEDHALOEVOLUTIONTREE_H_
#define DISTRIBUTEDHALOEVOLUTIONTREE_H_

#include "CosmologyToolsMacros.h"

// C++ STL includes
#include <vector>
#include <map>

// MPI
#include <mpi.h>

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
   * @brief get/set macro for MPI communicator.
   */
  GetNSetMacro(Communicator,MPI_Comm);

  /**
   * @brief Adds the given halos as nodes to this instance of the tree
   * @param halos the set of halos to be added in the tree.
   * @note Each tree node is identified by the halo hash-code.
   * @see Halo::GetHashCode()
   */
  void AppendNodes( std::vector< Halo >& halos );

  /**
   * @brief Creates an edge from the tree node corresponding to hashCode1
   * to the tree node corresponding to hashCode2.
   * @param hashCode1 the hashCode of the 1st tree node.
   * @param hashCode2 the hashCode of the 2nd tree node.
   * @pre this->HasNode( hashCode1 ) && this->HasNode( hashCode2 ) is true.
   */
  void CreateEdge(
          std::string hashCode1,
          std::string hashCode2 );

  /**
   * @brief Checks if the tree is empty.
   * @return true iff the tree is empty, else, false.
   */
  bool IsEmpty();

  /**
   * @brief Checks if the halo corresponding to the given hashCode exists
   * (in this process) on the tree.
   * @param hashCode the hash code corresponding to the halo in query.
   * @return true iff the halo exists in the tree, else false.
   */
  bool HasNode( std::string hashCode1 );

protected:
  std::map< std::string, Halo > Nodes; // List of nodes in the halo
  std::vector< std::string >    Edges; // List of edges between halos

  MPI_Comm Communicator;

  /**
   * @brief Relabels s.t. all nodes are uniquely identified by an integer ID.
   */
  void RelabelTreeNodes();

private:
  DISABLE_COPY_AND_ASSIGNMENT(DistributedHaloEvolutionTree);

};

} /* namespace cosmotk */
#endif /* DISTRIBUTEDHALOEVOLUTIONTREE_H_ */
