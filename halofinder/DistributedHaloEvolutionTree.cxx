#include "DistributedHaloEvolutionTree.h"

// CosmologyTools includes
#include "GenericIO.h"
#include "GenericIODefinitions.hpp"
#include "Halo.h"
#include "MergerTreeEvent.h"
#include "MergerTreeFileFormat.h"

// STL includes
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <set>
#include <limits>

//------------------------------------------------------------------------------

namespace cosmotk
{

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::DistributedHaloEvolutionTree()
{
  this->Communicator     = MPI_COMM_NULL;
  this->NumberOfMergers  = 0;
  this->NumberOfRebirths = 0;
  this->NumberOfSplits   = 0;
  this->NumberOfZombies  = 0;
  this->MergerTreeFileFormat = cosmotk::MergerTreeFileFormat::GENERIC_IO_POSIX;
}

//------------------------------------------------------------------------------
DistributedHaloEvolutionTree::~DistributedHaloEvolutionTree()
{
  this->Clear();
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::InsertNode(
      const HaloInfo &halo, unsigned char mask)
{
  assert("pre: corrupted merger-tree" && this->EnsureArraysAreConsistent());

  // STEP 0: Get hash code for node
  std::string hashCode = Halo::GetHashCodeForHalo(halo.Tag,halo.TimeStep);
  assert("pre: Encountered duplicate tree node!" && !this->HasNode(hashCode));

  // STEP 1: Insert node to list and create node-to-index mapping
  this->Nodes.push_back( halo );
  this->Node2Idx[ hashCode ] = this->Nodes.size()-1;

  // STEP 2: Initialize progenitor and descendant lists for the node
  std::vector< int > myProgenitors;
  myProgenitors.resize( 0 );
  std::vector< int > myDescendants;
  myDescendants.resize( 0 );
  this->Progenitors.push_back( myProgenitors );
  this->Descendants.push_back( myDescendants );

  // STEP 3: Store node event bitmask
  this->EventBitMask.push_back( mask );

  // STEP 4: Update various statistics
  this->UpdateNodeCounter(halo.TimeStep);

  if( MergerTreeEvent::IsEvent(mask,MergerTreeEvent::DEATH) )
    {
    ++this->NumberOfZombies;
    this->UpdateZombieCounter(halo.TimeStep);
    }

  if( MergerTreeEvent::IsEvent(mask,MergerTreeEvent::MERGE))
    {
    ++this->NumberOfMergers;
    }

  if( MergerTreeEvent::IsEvent(mask,MergerTreeEvent::REBIRTH) )
    {
    ++this->NumberOfRebirths;
    }

  if( MergerTreeEvent::IsEvent(mask,MergerTreeEvent::SPLIT) )
    {
    ++this->NumberOfSplits;
    }

  assert("post: corrupted merger-tree" && this->EnsureArraysAreConsistent());
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::LinkHalos(
      Halo* progenitor, Halo* descendant)
{
  // Sanity Checks
  assert("pre: corrupted merger-tree" && this->EnsureArraysAreConsistent());
  assert("pre: progenitor is NULL!" && (progenitor != NULL) );
  assert("pre: descendant is NULL!" && (descendant != NULL) );
  assert("pre: timesteps are out-of-synch!" &&
         (progenitor->TimeStep < descendant->TimeStep) );

  // STEP 0: If the progenitor is local, i.e., exists in this rank, update
  // its descendant list.
  if(this->HasNode(progenitor->GetHashCode()))
    {
    int idx = this->Node2Idx[ progenitor->GetHashCode() ];
    this->Descendants[ idx ].push_back( descendant->GlobalID );
    }

  // STEP 1: If the descendant is local, i.e., exists in this rank, update
  // its progenitor list.
  if(this->HasNode(descendant->GetHashCode()))
    {
    int idx = this->Node2Idx[ descendant->GetHashCode() ];
    this->Progenitors[ idx ].push_back( progenitor->GlobalID );
    }

  assert("post: corrupted merger-tree" && this->EnsureArraysAreConsistent());
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::Clear()
{
  this->Nodes.clear();
  this->Progenitors.clear();
  this->Descendants.clear();
  this->EventBitMask.clear();
  this->Node2Idx.clear();
  this->NodeCounter.clear();
  this->ZombieCounter.clear();
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfNodes(const int tstep)
{
  int N = 0;
  if(this->NodeCounter.find(tstep) != this->NodeCounter.end() )
    {
    N = this->NodeCounter[ tstep ];
    }
  return( N );
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetNumberOfZombieNodes(
    const int tstep)
{
  int N = 0;
  if( this->ZombieCounter.find(tstep) != this->ZombieCounter.end() )
    {
    N = this->ZombieCounter[ tstep ];
    }
  return( N );
}

//------------------------------------------------------------------------------
bool DistributedHaloEvolutionTree::HasNode(std::string hashCode)
{
  if( this->Node2Idx.find(hashCode) == this->Node2Idx.end() )
    {
    return false;
    }
  return true;
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::UpdateNodeCounter(const int tstep)
{
  if( this->NodeCounter.find(tstep) != this->NodeCounter.end() )
    {
    this->NodeCounter[tstep]++;
    }
  else
    {
    this->NodeCounter[tstep] = 1;
    }
}

//-----------------------------------------------------------------------------
void DistributedHaloEvolutionTree::UpdateZombieCounter(const int tstep)
{
  if( this->ZombieCounter.find(tstep) != this->ZombieCounter.end() )
    {
    this->ZombieCounter[tstep]++;
    }
  else
    {
    this->ZombieCounter[tstep] = 1;
    }
}

//------------------------------------------------------------------------------
int DistributedHaloEvolutionTree::GetTotalNumberOfBytes()
{
  // STEP 0: Compute local number of bytes
  int localNumBytes = 0;

  // tree nodes
  localNumBytes += this->Nodes.size()*sizeof(HaloInfo);

  // bit mask for each node
  localNumBytes += this->EventBitMask.size()*sizeof(unsigned char);

  // Sum progenitor list
  for(unsigned int idx=0; idx < this->Progenitors.size(); ++idx)
    {
    localNumBytes += this->Progenitors[idx].size()*sizeof(int);
    } // END for all progenitors

  // Sum descendant list
  for(unsigned int idx=0; idx < this->Descendants.size(); ++idx)
    {
    localNumBytes += this->Descendants[idx].size()*sizeof(int);
    }

  // Sum memory for Node2Idx map
  std::map< std::string, int >::iterator iter = this->Node2Idx.begin();
  for(; iter != this->Node2Idx.end(); ++iter )
    {
    localNumBytes += iter->first.size()*sizeof(char);
    localNumBytes += sizeof(int);
    } // END for for all items in the map

  // various statistics
  localNumBytes += 4*sizeof(int)+
      2*this->ZombieCounter.size()*sizeof(int)+
      2*this->NodeCounter.size()*sizeof(int);


  // STEP 1: Sum-reduce the local sum to compute the total across all ranks
  int totalNumBytes = 0;
  MPI_Allreduce(
      &localNumBytes,&totalNumBytes,1,MPI_INT,MPI_SUM,this->Communicator);
  return( totalNumBytes );
}

//------------------------------------------------------------------------------
void DistributedHaloEvolutionTree::WriteTree(std::string fileName)
{
  // STEP 0: Allocate writer
  GenericIO *writer = NULL;
  switch(this->MergerTreeFileFormat)
    {
    case MergerTreeFileFormat::GENERIC_IO_MPI:
      writer = new GenericIO(
          this->Communicator,fileName,GenericIO::FileIOMPI);
      break;
    case MergerTreeFileFormat::GENERIC_IO_POSIX:
      writer = new GenericIO(
          this->Communicator,fileName,GenericIO::FileIOPOSIX);
      break;
    default:
      std::cerr << "WARNING: invalid file format! ";
      std::cerr << "Defaulting to GENERIC_IO_POSIX\n";
      writer = new GenericIO(
          this->Communicator,fileName,GenericIO::FileIOPOSIX);
    } // END switch

  // STEP 1: Allocate temporary arrays for writting
  int N = this->GetNumberOfNodes();
  std::vector< ID_T > treeNodeIds(N+(CRCSize/sizeof(ID_T)),-1);
  std::vector< ID_T > haloTags(N+(CRCSize/sizeof(ID_T)),-1);
  std::vector< REAL > haloMass(N+(CRCSize/sizeof(REAL)),0.0);
  std::vector< int > tsteps(N+(CRCSize/sizeof(int)),-1);
  std::vector< REAL > redshift(N+(CRCSize/sizeof(REAL)),0.0);
  std::vector< POSVEL_T > center_x(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > center_y(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > center_z(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > mcx(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > mcy(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > mcz(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > vx(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > vy(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< POSVEL_T > vz(N+(CRCSize/sizeof(POSVEL_T)),0.0);
  std::vector< ID_T > descendant(N+(CRCSize/sizeof(ID_T)),-1);
  std::vector< unsigned char > eventMask(
      N+(CRCSize/sizeof(unsigned char)),0xFE /* poison */);

  // STEP 2: Fill arrays
  for(int i=0; i < N; ++i)
    {
    treeNodeIds[ i ] = this->Nodes[ i ].GlobalID;
    haloTags[ i ]    = this->Nodes[ i ].Tag;
    haloMass[ i ]    = this->Nodes[ i ].HaloMass;
    tsteps[ i ]      = this->Nodes[ i ].TimeStep;
    redshift[ i ]    = this->Nodes[ i ].Redshift;
    center_x[ i ]    = this->Nodes[ i ].Center[ 0 ];
    center_y[ i ]    = this->Nodes[ i ].Center[ 1 ];
    center_z[ i ]    = this->Nodes[ i ].Center[ 2 ];
    mcx[ i ]         = this->Nodes[ i ].MeanCenter[ 0 ];
    mcy[ i ]         = this->Nodes[ i ].MeanCenter[ 1 ];
    mcz[ i ]         = this->Nodes[ i ].MeanCenter[ 2 ];
    vx[ i ]          = this->Nodes[ i ].AverageVelocity[ 0 ];
    vy[ i ]          = this->Nodes[ i ].AverageVelocity[ 1 ];
    vz[ i ]          = this->Nodes[ i ].AverageVelocity[ 2 ];
    descendant[ i ]  = this->GetDescendant( i );
    eventMask[ i ]   = this->EventBitMask[ i ];
    } // END for all nodes

  // STEP 3: Register variables & data arrays to the writer
  writer->addVariable(
      "tree_node_id",&treeNodeIds[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("halo_tag",&haloTags[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("halo_mass",&haloMass[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("timestep",&tsteps[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("redshift",&redshift[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("center_x",&center_x[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("center_y",&center_y[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("center_z",&center_z[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("mean_center_x",&mcx[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("mean_center_y",&mcy[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("mean_center_z",&mcz[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("v_x",&vx[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("v_y",&vy[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("v_z",&vz[0],GenericIO::VarHasExtraSpace);
  writer->addVariable(
      "descendant_id",&descendant[0],GenericIO::VarHasExtraSpace);
  writer->addVariable("event_mask",&eventMask[0],GenericIO::VarHasExtraSpace);

  // STEP 4: Write the data
  writer->write();
  delete writer;

  // STEP 5: De-allocate temporary arrays & data-structures
  treeNodeIds.clear();
  haloTags.clear();
  haloMass.clear();
  tsteps.clear();
  redshift.clear();
  center_x.clear();
  center_y.clear();
  center_z.clear();
  mcx.clear();
  mcy.clear();
  mcz.clear();
  vx.clear();
  vy.clear();
  vz.clear();
  descendant.clear();
  eventMask.clear();
}

//------------------------------------------------------------------------------
ID_T DistributedHaloEvolutionTree::GetDescendant(const int i)
{
  assert("pre: local index out-of-bounds!" &&
          (i >= 0) && (i < this->Nodes.size() ) );

  if( this->Descendants[ i ].size() == 0 )
    {
    return( -1 );
    }
  else if( this->Descendants[ i ].size() > 1 )
    {
    // TODO: must select descendant with highest mass
    return( this->Descendants[ i ][ 0 ]);
    }
  assert("pre: No single descendant!" && (this->Descendants[i].size()==1));
  return( this->Descendants[ i ][ 0 ]);
}

//------------------------------------------------------------------------------
std::string DistributedHaloEvolutionTree::ToString()
{
  assert("pre: corrupted merger-tree" && this->EnsureArraysAreConsistent());

  std::ostringstream oss;
  oss << "NUMNODES " << this->GetNumberOfNodes() << std::endl;
  oss << "ID\t";
  oss << "TIMESTEP\t";
  oss << "TAG\t";
  oss << "REDSHIFT\t";
  oss << "CENTER_X\t";
  oss << "CENTER_Y\t";
  oss << "CENTER_Z\t";
  oss << "MEAN_CENTER_X\t";
  oss << "MEAN_CENTER_Y\t";
  oss << "MEAN_CENTER_Z\t";
  oss << "V_X\t";
  oss << "V_Y\t";
  oss << "V_Z\t";
  oss << "MASS\t";
  oss << "DESCENDANTS\t";
  oss << "PROGENITORS\t";
  oss << "EVENT_TYPE\t";
  oss << std::endl;

  for(int i=0; i < this->GetNumberOfNodes(); ++i)
    {
    oss << std::scientific
        << std::setprecision(std::numeric_limits<POSVEL_T>::digits10);
    oss << this->Nodes[ i ].GlobalID             << "\t";
    oss << this->Nodes[ i ].TimeStep             << "\t";
    oss << this->Nodes[ i ].Tag                  << "\t";
    oss << this->Nodes[ i ].Redshift             << "\t";
    oss << this->Nodes[ i ].Center[0]            << "\t";
    oss << this->Nodes[ i ].Center[1]            << "\t";
    oss << this->Nodes[ i ].Center[2]            << "\t";
    oss << this->Nodes[ i ].MeanCenter[ 0 ]      << "\t";
    oss << this->Nodes[ i ].MeanCenter[ 1 ]      << "\t";
    oss << this->Nodes[ i ].MeanCenter[ 2 ]      << "\t";
    oss << this->Nodes[ i ].AverageVelocity[ 0 ] << "\t";
    oss << this->Nodes[ i ].AverageVelocity[ 1 ] << "\t";
    oss << this->Nodes[ i ].AverageVelocity[ 2 ] << "\t";
    oss << this->Nodes[ i ].HaloMass << "\t";

    oss << "[ ";
    for(int didx=0; didx < this->Descendants[i].size(); ++didx)
      {
      oss << this->Descendants[i][didx] << " ";
      } // END for all descendants

    // Put a -1 if there are no descendants. Note, this
    // should only happen at the very final time-step.
    if( this->Descendants[i].size() == 0 )
      {
      oss << "-1 ";
      }
    oss << "]\t";

    oss << "[ ";
    for(int pidx=0; pidx < this->Progenitors[i].size(); ++pidx)
      {
      oss << this->Progenitors[i][pidx] << " ";
      }

    // Put -999 here to indicate an empty progenitor list.
    if( this->Progenitors[i].size() == 0 )
      {
      oss << "-999 ";
      }
    oss << "]\t";

    unsigned char bitmask = this->EventBitMask[ i ];
    oss << MergerTreeEvent::GetEventString(bitmask);

    if( this->Descendants[i].size()==0 )
      {
      oss << "(**FINAL TIMESTEP**)";
      }

    oss << std::endl;
    } // END for all nodes

  return( oss.str() );
}


} /* namespace cosmotk */
