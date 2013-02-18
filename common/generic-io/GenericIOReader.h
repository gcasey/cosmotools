/**
 * @brief GenericIOReader defines the interface of the GenericIO readers.
 * This is a pure abstract class -- all methods are defined by concrete
 * classes.
 */
#ifndef GENERICIOREADER_H_
#define GENERICIOREADER_H_

#include "GenericIOBase.h" // Base class

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

private:
  DISABLE_COPY_AND_ASSIGNMENT(GenericIOReader);
};

} /* namespace cosmotk */
#endif /* GENERICIOREADER_H_ */
