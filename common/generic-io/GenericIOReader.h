/**
 * @brief GenericIOReader defines the interface of the GenericIO readers.
 * This is a pure abstract class -- all methods are defined by concrete
 * classes.
 */
#ifndef GENERICIOREADER_H_
#define GENERICIOREADER_H_

#include "GenericIOBase.h" // Base class
#include "GenericIODefinitions.hpp"

#include <map>    // for STL map
#include <string> // for STL string
#include <vector> // for STL vector

namespace cosmotk
{

class GenericIOReader : public GenericIOBase
{
public:
  GenericIOReader();
  virtual ~GenericIOReader();

  /**
   * @brief Return the variable size of the ith variable
   * @param i the index of the variable in query
   * @return varsize the size of the variable, e.g., 8 if it's a double, etc.
   */
  uint64_t GetVariableSize(const int i)
    {
    assert("pre: variable index is out-of-bounds!" &&
           (i >= 0) && (i < this->VH.size() ) );
    return( this->VH[ i ].Size );
    };

  /**
   * @brief Return the name of the ith variable
   * @param i the index of the variable in query
   * @return varname the variable name
   */
  std::string GetVariableName(const int i)
    {
    assert("pre: variable index is out-of-bounds!" &&
            (i >= 0) && (i < this->VH.size() ) );
    return( std::string( this->VH[ i ].Name ) );
    };

  /**
   * @brief Return the number of variables in the file
   * @return nvar the number of variables in the file
   */
  int GetNumberOfVariables()
    {return this->VH.size();};

  /**
   * @brief Return the number of assigned blocks
   * @return nassigned the number of assigned blocks to this process.
   */
  int GetNumberOfAssignedBlocks()
    {return this->AssignedBlocks.size();};

  /**
   * @brief Return the total number of blocks in the file.
   * @return nblocks the total number of blocks in the file.
   * @note The number of blocks in the file is equivalent to
   * the number of processes that wrote the file.
   */
  int GetTotalNumberOfBlocks()
    {return this->GH.NRanks;};

  /**
   * @brief Return the total number of elements.
   * @return N the total number of elements.
   */
  int GetTotalNumberOfElements()
    {return this->GH.NElems;};

  /**
   * @brief Checks if the given variable exists.
   * @param varName the variable name.
   * @return status true if the variable exists, else false.
   */
  bool HasVariable(std::string varName)
    { return( (this->GetVariableIndex(varName)!=-1) ); };

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
  bool SplitMode;

  GlobalHeader GH;
  std::vector< VariableHeader > VH;
  std::vector< RankHeader > RH;

  std::map<std::string, int> VariableName2IndexMap;

  // Stores the mapping of data-blocks to a file ID. Note, this is only
  // applicable iff this->SplitMode==true.
  std::map<int,int> BlockToFileMap;

  // The number of files the data has been split to. Note, this is only
  // applicable iff this->SplitMode==true.
  int NumberOfFiles;

  // List of blocks that will be read by this process. Recall, the number of
  // blocks in the file may not always match the number of ranks that this
  // reader instance is running. Hence, each rank may be assigned more than
  // one block in the file. The vector of AssignedBlocks holds the list of
  // blocks for this process.
  std::vector< int > AssignedBlocks;

  /**
   * @brief Checks if the underlying files has been written in split mode
   * @return status true, if the mode is split, else false.
   */
  bool IsSplitMode()
    { return( this->SplitMode ); };

  /**
   * @brief Based on the global & variable headers, this method determines if
   * the file has been written in split mode or not.
   * @pre The variables must have been indexed prior to calling this method,
   * i.e., this method should be called after IndexVariables() has been called.
   */
  void DetermineFileType();

  /**
   * @brief Builds an index based on variable name.
   */
  void IndexVariables();

  /**
   * @brief Returns the number of elements at the given block.
   * @param localBlkIdx the local block index of the block in query.
   * @return N the number of elements stored in the requested block.
   * @note Assumes that the rank header information has been read.
   */
  int GetNumberOfElementsForBlock(const int localBlkIdx);

  /**
   * @brief Return the index of the variable with the given name.
   * @param name the name of the variable in query
   * @return idx the corresponding index of the variable
   * @post idx >= 0, idx == -1 iff a variable is not found.
   */
  int GetVariableIndex(const std::string name);

  /**
   * @brief Computes the variable offset within the given global block idx.
   * @param vidx the variable index.
   * @param localBlkIdx the local block index.
   * @return offSet the variable offset.
   * @note This method uses the information in the rank header
   */
  uint64_t GetVariableOffSet(int vidx, int localBlkIdx);

private:
  DISABLE_COPY_AND_ASSIGNMENT(GenericIOReader);
};

} /* namespace cosmotk */
#endif /* GENERICIOREADER_H_ */
