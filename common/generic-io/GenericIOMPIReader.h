/**
 * @brief A concrete instance of GenericIOReader that implements the GenericIO
 * reader interface using MPI.
 */
#ifndef GENERICIOMPIREADER_H_
#define GENERICIOMPIREADER_H_

// CosmologyTools includes
#include "CosmologyToolsMacros.h"
#include "GenericIOReader.h"

// C/C++ includes
#include <sys/types.h> // for off_t


// MPI includes
#include <mpi.h>

namespace cosmotk
{

class GenericIOMPIReader : public GenericIOReader
{
public:
  GenericIOMPIReader();
  virtual ~GenericIOMPIReader();

  /**
   * @brief Opens and reads the header of the file and initializes internal
   * data-structures.
   * @pre !This->FileName.empty()
   * @pre this->Communicator != MPI_COMM_NULL
   */
  virtual void OpenAndReadHeader( bool skipBlockHeaders=false );

  /**
   * @brief Read the headers of each assigned block.
   * @see ReadBlockHeader
   */
  void ReadBlockHeaders();

  /**
   * @brief Overrides GetNumberOfElements to provide support for SplitMode
   * @return N the number of elements to read.
   * @see GenericIOReader::GetNumberOfElements()
   */
  virtual int GetNumberOfElements();

  /**
   * @brief Reads the data in to the user-supplied registered arrays.
   * @note The user should have registered the arrays to read via calls to
   * the AddVariable method.
   * @see GenericIOBase::AddVariable.
   * @pre this->Communicator != MPI_COMM_NULL
   */
  virtual void ReadData();

  /**
   * @brief Closes the file
   */
  virtual void Close();

protected:
  MPI_File FH;

  // Flag that indicates that all data requests should be proxied to the
  // internal readers.
  bool ProxyEnabled;

  /**
   * @brief List of internal readers
   */
  GenericIOMPIReader** InternalReaders;

  /**
   * @brief Reads in the BlockToFile mapping.
   * @pre this->SplitMode==true.
   */
  void ReadBlockToFileMap();

  /**
   * @brief Reads the header of the ith variable
   * @param vh the variable header where the dara will be read in.
   * @pre vh != NULL
   * @see ReadVariableHeaders
   */
  void ReadVariableHeader( const int i, VariableHeader& vh );

  /**
   * @brief Reads in the variable headers for the number of variables.
   * @pre This method assumes that the global header has been read in.
   * @see ReadVariableHeader
   */
  void ReadVariableHeaders();

  /**
   * @brief Reads the global and variable header of the file and broadcasts
   * them to all ranks.
   */
  void ReadHeader();

  /**
   * @brief Reads the block header corresponding to the given block index.
   * @param blkIdx the index of the block to reader
   * @param blockHeader data-structure where to read in the block header
   * @see ReadBlockHeaders
   */
  void ReadBlockHeader(const int blkIdx, RankHeader& blockHeader);

  /**
   * @brief If SplitMode is used, this method will attach a reader for each
   * corresponding file.
   */
  void SetupInternalReaders();

  /**
   * @brief Read the data in split mode.
   */
  void ReadSplitModeData();

  /**
   * @brief Reads data from a single file.
   */
  void ReadSingleFileData();

  /**
   * @brief Reads data into the user-supplied buffer from the MPI file handle
   * @param buf the buffer where the data will be read into
   * @param count the number of bytes to read
   * @param offset the offset from which to read
   * @param variableName name of the data being read, primarily, used for
   * debugging and error reporting.
   */
  virtual void Read(void *buf, size_t count,
                      off_t offset, const std::string &variableName );

private:
  DISABLE_COPY_AND_ASSIGNMENT(GenericIOMPIReader);
};

} /* namespace cosmotk */
#endif /* GENERICIOMPIREADER_H_ */
