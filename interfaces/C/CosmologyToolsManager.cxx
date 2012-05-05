#include "CosmologyToolsManager.h"


namespace cosmologytools
{

CosmologyToolsManager::CosmologyToolsManager()
{
  this->Layout                = MemoryLayout::ROWMAJOR;
  this->Communicator          = MPI_COMM_WORLD;
  this->EnableVis             = false;
  this->ArrayStride           = 6;
  this->HaloTrackingFrequency = 5;
  this->HaloFinder            = 0;
}

//-----------------------------------------------------------------------------
CosmologyToolsManager::~CosmologyToolsManager()
{

}

} /* namespace cosmologytools */
