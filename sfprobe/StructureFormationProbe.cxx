#include "StructureFormationProbe.h"

// C++ includes
#include <cstdlib>
#include <cassert>
#include <fstream>
#include <sstream>

#define IMIN(ext) ext[0]
#define IMAX(ext) ext[1]
#define JMIN(ext) ext[2]
#define JMAX(ext) ext[3]
#define KMIN(ext) ext[4]
#define KMAX(ext) ext[5]

namespace cosmologytools {

StructureFormationProbe::StructureFormationProbe()
{
  this->Initialize();
}

//------------------------------------------------------------------------------
StructureFormationProbe::~StructureFormationProbe()
{
  if( this->Connectivity != NULL )
   {
   delete [] this->Connectivity;
   }
}

//------------------------------------------------------------------------------
void StructureFormationProbe::Initialize()
{
  this->Connectivity = NULL;
  this->Particles    = NULL;
  this->NumTets      = 0;
  for( int i=0; i < 3; ++i )
    {
    this->WholeExtent[i*2]   = this->Extent[ i*2 ]   = 0;
    this->WholeExtent[i*2+1] = this->Extent[ i*2+1 ] = 4;
    }
}

//------------------------------------------------------------------------------
void StructureFormationProbe::SetGridExtent(int extent[6])
{
  assert( "pre: grid extent must be 3-D" && this->Is3DExtent(extent) );

  for( int i=0; i < 6; ++i )
    {
    this->Extent[i] = extent[i];
    }

  if( this->Connectivity != NULL )
    {
    delete [] this->Connectivity;
    }

  int numVoxels = this->GetNumberOfCells( extent );
  this->NumTets   = 5*numVoxels;
  this->Connectivity = new int[ 4*this->NumTets ];
}

//------------------------------------------------------------------------------
void StructureFormationProbe::Tesselate()
{
  this->TesselateGridExtent();
}

#define ADDTET(tetIdx, V0, V1, V2, V3 )   \
    this->Connectivity[tetIdx*4]   = V0;  \
    this->Connectivity[tetIdx*4+1 ]= V1;  \
    this->Connectivity[tetIdx*4+2 ]= V2;  \
    this->Connectivity[tetIdx*4+3 ]= V3;  \
    ++tetIdx;

//------------------------------------------------------------------------------
void StructureFormationProbe::TesselateGridExtent()
{
  assert("pre: tetrahedral connectivity array is NULL" &&
         (this->Connectivity != NULL) );

  int V[8];
  int tetIdx  = 0;
  int cellIdx = 0;
  for( int i=IMIN(this->Extent); i < IMAX(this->Extent); ++i )
    {
    for( int j=JMIN(this->Extent); j < JMAX(this->Extent); ++j )
      {
      for( int k=KMIN(this->Extent); k < KMAX(this->Extent); ++k )
        {
        cellIdx = this->GetGlobalLinearIndex(i,j,k);

        // Nodes on voxel base
        V[0] = this->GetLinearIndex(i,j,k);
        V[1] = this->GetLinearIndex(i+1,j,k);
        V[2] = this->GetLinearIndex(i+1,j+1,k);
        V[3] = this->GetLinearIndex(i,j+1,k);

        // Nodes on voxel top
        V[4] = this->GetLinearIndex(i,j,k+1);
        V[5] = this->GetLinearIndex(i+1,j,k+1);
        V[6] = this->GetLinearIndex(i+1,j+1,k+1);
        V[7] = this->GetLinearIndex(i,j+1,k+1);

        // Set tetrahedral connectivity
        if( cellIdx % 2 )
          {
          ADDTET(tetIdx,V[1],V[2],V[3],V[6]);
          ADDTET(tetIdx,V[0],V[1],V[3],V[4]);
          ADDTET(tetIdx,V[4],V[6],V[5],V[1]);
          ADDTET(tetIdx,V[4],V[6],V[7],V[3]);
          ADDTET(tetIdx,V[3],V[6],V[4],V[1]);
          }
        else
          {
          ADDTET(tetIdx,V[1],V[5],V[2],V[0]);
          ADDTET(tetIdx,V[2],V[3],V[0],V[7]);
          ADDTET(tetIdx,V[2],V[5],V[6],V[7]);
          ADDTET(tetIdx,V[0],V[7],V[4],V[5]);
          ADDTET(tetIdx,V[0],V[2],V[7],V[5]);
          }

        } // END for all k
      } // END for all j
    } // END for all i
}

//------------------------------------------------------------------------------
void StructureFormationProbe::WriteTesselation(char* fileName)
{
  std::ofstream ofs;
  ofs.open(fileName);
  assert("pre: Cannot open file!" && ofs.is_open() );

  // Write VTK header
  ofs << "# vtk DataFile Version 3.0\n";
  ofs << "Tesselation\n";
  ofs << "ASCII\n";
  ofs << "DATASET UNSTRUCTURED_GRID\n";

  // Write points
  ofs << "POINTS " << this->NumParticles << " double\n";
  for( int i=0; i < this->NumParticles; ++i )
    {
    ofs << this->Particles[i*this->Stride]   << " ";
    ofs << this->Particles[i*this->Stride+1] << " ";
    ofs << this->Particles[i*this->Stride+2] << std::endl;
    }

  // Write tets
  std::ostringstream constream;
  std::ostringstream typestream;
  std::ostringstream idxstream;
  for( int i=0; i < this->NumTets; ++i )
    {
    constream  << "4 " << this->Connectivity[i*4] << " ";
    constream  << this->Connectivity[i*4+1] << " ";
    constream  << this->Connectivity[i*4+2] << " ";
    constream  << this->Connectivity[i*4+3] << std::endl;
    typestream << "10\n";
    idxstream  << i << std::endl;
    } // END for all tets

  ofs << "CELLS " << this->NumTets << " " << this->NumTets*5 << std::endl;
  ofs << constream.str();

  ofs << "CELL_TYPES " << this->NumTets << std::endl;
  ofs << typestream.str();

  ofs << "CELL_DATA " << this->NumTets << std::endl;
  ofs << "SCALARS CellIdx int 1\n";
  ofs << "LOOKUP_TABLE default\n";
  ofs << idxstream.str();
  ofs.close();
}

//------------------------------------------------------------------------------
int StructureFormationProbe::GetGlobalLinearIndex(
    const int i, const int j,const int k)
{
  // Get Extent dimensions
  int N2  = JMAX(this->WholeExtent)-JMIN(this->WholeExtent)+1;
  int N1  = IMAX(this->WholeExtent)-IMIN(this->WholeExtent)+1;

  // Return the linear index
  return( (k*N2+j)*N1+i );
}

//------------------------------------------------------------------------------
int StructureFormationProbe::GetLinearIndex(const int i,const int j,const int k)
{
  // Get local indices staring from 0
  int li = i-IMIN(this->Extent);
  int lj = j-JMIN(this->Extent);
  int lk = k-KMIN(this->Extent);

  // Get Extent dimensions
  int N2  = JMAX(this->Extent)-JMIN(this->Extent)+1;
  int N1  = IMAX(this->Extent)-IMIN(this->Extent)+1;

  // Return linear index
  return( (lk*N2+lj)*N1+li );
}

//------------------------------------------------------------------------------
bool StructureFormationProbe::Is3DExtent( int ext[6] )
{
  bool status = false;
  if( (IMAX(ext) > IMIN(ext) ) &&
      (JMAX(ext) > JMIN(ext) ) &&
      (KMAX(ext) > KMIN(ext)) )
    {
    status = true;
    }
  return( status );
}

//------------------------------------------------------------------------------
int StructureFormationProbe::GetNumberOfNodes( int ext[6] )
{
  int N = 0;
  int idim = IMAX(ext)-IMIN(ext)+1;
  int jdim = JMAX(ext)-JMIN(ext)+1;
  int kdim = KMAX(ext)-KMIN(ext)+1;
  N        = idim * jdim * kdim;
  return( N );
}

//------------------------------------------------------------------------------
int StructureFormationProbe::GetNumberOfCells( int ext[6] )
{
  int N = 0;
  int idim = IMAX(ext)-IMIN(ext);
  int jdim = JMAX(ext)-JMIN(ext);
  int kdim = KMAX(ext)-KMIN(ext);
  N        = idim * jdim * kdim;
  return( N );
}

} /* namespace hacctools */
