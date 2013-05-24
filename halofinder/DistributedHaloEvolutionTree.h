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

// Forward declarations
struct HaloInfo;


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

  // Get macros for some of the statistics
  GetMacro(NumberOfMergers,int);
  GetMacro(NumberOfRebirths,int);
  GetMacro(NumberOfSplits,int);
  GetMacro(NumberOfZombies,int);

  /**
   * @brief Inserts the given halo as a node in the tree.
   * @param halo the halo to add
   */
  void InsertNode(const HaloInfo &halo, unsigned char mask);

  /**
   * @brief Links the progenitor from a previous timestep to the descendant
   * halo in the current timestep.
   * @param progenitor pointer to the progenitor.
   * @param descendant pointer to the descendant.
   * @pre progenitor != NULL.
   * @pre descendant != NULL.
   */
  void LinkHalos(Halo* progenitor, Halo* descendant);

  /**
   * @brief Checks if the tree is empty.
   * @return true iff the tree is empty, else, false.
   */
  bool IsEmpty()
    { return( this->Nodes.empty() ); };

  /**
    * @brief Returns the total number of nodes in the tree on this process.
    * @return N the local number of nodes.
    */
   int GetNumberOfNodes()
     {return this->Nodes.size();};

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
   * @brief Checks if the halo corresponding to the given hashCode exists
   * (in this process) on the tree.
   * @param hashCode the hash code corresponding to the halo in query.
   * @return true iff the halo exists in the tree, else false.
   */
  bool HasNode( std::string hashCode );

  /**
   * @brief Clears out all the data in this tree
   * @post this->IsEmpty()==true
   */
  void Clear();

  /**
   * @brief Returns the number of bytes used by this instance, on all ranks.
   * @return N the total number of bytes.
   * @note This method is collective, hence, it must be called by all ranks.
   */
  int GetTotalNumberOfBytes();

  /**
   * @brief Barrier synchronization among all processes
   */
  void Barrier() { MPI_Barrier(this->Communicator); }

  /**
   * @brief Creates an ASCII string representation of the given tree
   * @return a string representing the data in the tree
   * @note Used mainly for debugging.
   */
  std::string ToString();

  /**
   * @brief Writes the Tree in GenericIO.
   */
  void WriteTree( std::string fileName );

protected:

  // User-supplied MPI communicator
  MPI_Comm Communicator;

  // List of tree nodes
  std::vector<HaloInfo> Nodes; // Stores information for each tree node

  // For each node, store a list of progenitors
  std::vector< std::vector<int> > Progenitors;

  // For each node, store a list of descendants
  std::vector< std::vector<int> > Descendants;

  // For each node, store the bitmask that encodes the various events.
  std::vector<unsigned char> EventBitMask;

  // Mapping of halo hashcode to the index of the halo in the Nodes vector.
  // Note, Halo hash codes are generated using Halo::GetHashCodeForHalo method.
  // This data-structure is used to determine if a node is in the tree.
  std::map< std::string, int > Node2Idx;

  // VARIOUS STATISTICS
  // Number of Nodes at each time-step
  std::map<int,int> NodeCounter;

  // Number of Zombies at each time-step
  std::map<int,int> ZombieCounter; // Counts number of zombies at a given tstep


  // Total counts
  int NumberOfZombies;
  int NumberOfMergers;
  int NumberOfRebirths;
  int NumberOfSplits;

  /**
   * @brief Returns the descendant of the ith node.
   * @param i the local node index, i.e., index to the Nodes array.
   * @return idx the ID of the descendant.
   * @note Most times a node will have a single descendant. However, when
   * SPLIT events occur, it is possible to have more than one descendant.
   * In this case, the method will select a descendant with the highest mass.
   */
  ID_T GetDescendant(const int i);

  /**
   * @brief Ensures that the data arrays used internally are consistent.
   * @return status true if consistent, else, false.
   * @note Used mostly for sanity checking.
   */
  bool EnsureArraysAreConsistent()
    {
    int N = this->GetNumberOfNodes();
    return( (N==this->Progenitors.size()) &&
             (N==this->Descendants.size()) &&
             (N==this->EventBitMask.size()));
    }

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

private:
  DISABLE_COPY_AND_ASSIGNMENT(DistributedHaloEvolutionTree);

};

} /* namespace cosmotk */
#endif /* DISTRIBUTEDHALOEVOLUTIONTREE_H_ */
