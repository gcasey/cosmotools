/**
 * @brief GenericIOReader defines the interface of the GenericIO readers.
 * This is a pure abstract class -- all methods are defined by concrete
 * classes.
 */
#ifndef GENERICIOREADER_H_
#define GENERICIOREADER_H_

#include "GenericIOBase.h" // Base class
#include "GenericIODefinitions.hpp"

#include <vector> // for STL vector
#include <map>    // for STL map

namespace cosmotk
{

class GenericIOReader : public GenericIOBase
{
public:
  GenericIOReader();
  virtual ~GenericIOReader();

  /**
   * @brief Opens and reads the header of the file and initializes internal
   * data-structures.
   * @pre !This->FileName.empty()
   */
  virtual void OpenAndReadHeader()  = 0;

  /**
   * @brief Returns the number of elements that will be read (by this process)
   * @note All variables have the same size.
   * @return N the number of elements to read
   * @post N >= 0.
   */
  virtual int GetNumberOfElements() = 0;

  /**
   * @brief Reads the data in to the user-supplied registered arrays.
   * @note The user should have registered the arrays to read via calls to
   * the AddVariable method.
   * @see GenericIOBase::AddVariable.
   */
  virtual void ReadData() = 0;

  /**
   * @brief Closes the file
   */
  virtual void Close() = 0;

protected:
  bool SwapEndian;

  GlobalHeader GH;
  std::vector< VariableHeader > VH;
  std::vector< RankHeader > RH;

  std::map<std::string, int> VariableName2IndexMap;

  // List of blocks that will be read by this process. Recall, the number of
  // blocks in the file may not always match the number of ranks that this
  // reader instance is running. Hence, each rank may be assigned more than
  // one block in the file. The vector of AssignedBlocks holds the list of
  // blocks for this process.
  std::vector< int > AssignedBlocks;

  /**
   * @brief Builds an index based on variable name.
   */
  void IndexVariables();

private:
  DISABLE_COPY_AND_ASSIGNMENT(GenericIOReader);
};

} /* namespace cosmotk */
#endif /* GENERICIOREADER_H_ */
