#include "GenericIOReader.h"

#include <iostream>

namespace cosmotk
{

//------------------------------------------------------------------------------
GenericIOReader::GenericIOReader()
{
  this->SwapEndian=false;
}

//------------------------------------------------------------------------------
GenericIOReader::~GenericIOReader()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void GenericIOReader::IndexVariables()
{
  std::string name;
  for(unsigned int i=0; i < this->VH.size(); ++i)
    {
    name = std::string( this->VH[i].Name );
    this->VariableName2IndexMap[ name ]= i;
    }
}

//------------------------------------------------------------------------------
int GenericIOReader::GetVariableIndex(const std::string name)
{
  int idx = -1;
  if( this->VariableName2IndexMap.find(name) !=
       this->VariableName2IndexMap.end())
    {
    idx = this->VariableName2IndexMap[name];
    }
  else
    {
    std::cerr << "WARNING: No index for variable " << name << "!\n";
    }
  return( idx );
}

//------------------------------------------------------------------------------
int GenericIOReader::GetNumberOfElementsForBlock(const int blkidx)
{
  assert("pre: blkIdx is out-of-bounds!" &&
          (blkidx >= 0) && (blkidx < this->AssignedBlocks.size() ) );
  assert("pre: blkIdx is out-of-bounds!" &&
          (blkidx >= 0) && (blkidx < this->RH.size() ) );
  return( this->RH[blkidx].NElems);
}

//------------------------------------------------------------------------------
uint64_t GenericIOReader::GetVariableOffSet(int vidx, int localBlkIdx)
{
  // Sanity checks!
  assert("pre: variable index is out-of-bounds!" &&
       (vidx >= 0) && (vidx < this->GH.NVars) );
  assert("pre: local block index is out-of-bounds!" &&
       (localBlkIdx >= 0) && (localBlkIdx < this->RH.size() ) );


  int N = this->GetNumberOfElementsForBlock(localBlkIdx);
  uint64_t offSet = this->RH[ localBlkIdx ].Start;
  for(int i=0; i < vidx; ++i)
    {
    offSet += N*this->VH[i].Size + CRCSize;
    }
  return( offSet );
}

} /* namespace cosmotk */
