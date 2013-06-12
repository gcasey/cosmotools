#include "ParallelStructureFormationProbe.h"
#include "ExtentPartitioner.h"

#include <sstream>

namespace cosmologytools {

ParallelStructureFormationProbe::ParallelStructureFormationProbe()
{
//  this->Initialize();
  this->Communicator = MPI_COMM_WORLD;
}

//------------------------------------------------------------------------------
ParallelStructureFormationProbe::~ParallelStructureFormationProbe()
{
  // TODO Auto-generated destructor stub
}

//------------------------------------------------------------------------------
void ParallelStructureFormationProbe::Tesselate()
{
//  // STEP 0: Get number of ranks, N, and the rank of this process
//  int N      = this->GetNumberOfRanks();
//  int myRank = this->Rank();
//
//  // STEP 1: Partition the Langrangian (initial) space to processes
//  ExtentPartitioner partitioner;
//  partitioner.SetGlobalExtent( this->WholeExtent );
//  partitioner.SetNumberOfPartitions( N );
//  partitioner.Partition();
//  this->Barrier();
//
//  int ext[6];
//  partitioner.GetExtent( myRank, ext );
//  this->SetGridExtent( ext );
//
//  // STEP 2: Tesselate
//  this->TesselateGridExtent();
//
//  std::ostringstream oss;
//  oss << "InitialTriangulation_" << myRank << ".vtk";
//  this->WriteTesselation( const_cast<char*>(oss.str().c_str()) );
//  this->Barrier();
}

} /* namespace hacctools */
