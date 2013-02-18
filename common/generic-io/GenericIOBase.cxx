#include "GenericIOBase.h"

namespace cosmotk
{

GenericIOBase::GenericIOBase()
{
  this->FileName   = "";
  this->Comm       = MPI_COMM_NULL;
  this->IOStrategy = FileIOUndefined;
}

//------------------------------------------------------------------------------
GenericIOBase::~GenericIOBase()
{
  // TODO Auto-generated destructor stub
}

} /* namespace cosmotk */
