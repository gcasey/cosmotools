#include "GenericIOBase.h"

namespace cosmotk
{

GenericIOBase::GenericIOBase()
{
  this->FileName   = "";
  this->IOStrategy = FileIOUndefined;
}

//------------------------------------------------------------------------------
GenericIOBase::~GenericIOBase()
{
  // TODO Auto-generated destructor stub
}

} /* namespace cosmotk */
