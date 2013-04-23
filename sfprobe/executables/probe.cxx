/**
 * A simple program to test the functionality of the structure formation probe
 */

// C/C++ includes
#include <cassert>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>

// CosmologyTools macros
#include "CosmologyToolsMacros.h"
#include "StructureFormationProbe.h"

// VTK includes
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCellTypes.h"
#include "vtkDoubleArray.h"
#include "vtkIntArray.h"
#include "vtkPCosmoReader.h"
#include "vtkPointData.h"
#include "vtkUniformGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"
#include "vtkXMLImageDataWriter.h"

#ifdef USEDAX
 // Dax includes
#  define DAX_DEVICE_ADAPTER DAX_DEVICE_ADAPTER_TBB
#  include <dax/CellTag.h>
#  include <dax/cont/arg/ExecutionObject.h>
#  include <dax/cont/DeviceAdapter.h>
#  include <dax/cont/Scheduler.h>
#  include <dax/cont/Timer.h>
#  include <dax/Extent.h>
#  include "dax/MapTetsToPoints.h"
#  include "dax/TetStructure.h"

#endif


//------------------------------------------------------------------------------
//  GLOBAL DEFINITIONS
//------------------------------------------------------------------------------
int NSTEPS = 1;
std::string files[]= {
 "initial.z20.cosmo"
};

//std::string files[]= {
// "initial.z20.cosmo",
// "initial.z10.cosmo",
// "initial.z5.cosmo",
// "initial.z2.cosmo",
// "initial.z0.cosmo"
//};

/**
 * @brief A simple object to store program parameters
 */
struct ProgramParams {
  std::string BaseDir; // base directory where the
  REAL rL;             // The domain length
  INTEGER LDIM;        // dimensions for the langrangian grid
  INTEGER GDIM;        // dimensions for the probe grid
  INTEGER Fringe;      // Fringe parameter used to guard from periodic bndries

  void Println()
    {
    std::cout << "Base directory: " << this->BaseDir << std::endl;
    std::cout << "rL: " << this->rL << std::endl;
    std::cout << "LDIM: " << this->LDIM << std::endl;
    std::cout << "GDIM: " << this->GDIM << std::endl;
    std::cout << "FRINGE: " << this->Fringe << std::endl;
    std::cout.flush();
    }
} Parameters;


/**
 * @class Particles
 * @brief A simple class to store the particle positions.
 */
class SimParticles
{
public:
  REAL *Positions;
  INTEGER *GlobalIDs;
  INTEGER N;

  SimParticles()
  {
    this->GlobalIDs = NULL;
    this->Positions = NULL;
    this->N = 0;
  }

  ~SimParticles()
  {
    this->Delete();
  }

  void Allocate(const INTEGER n)
  {
    this->N = n;
    this->Positions = new REAL[3*this->N];
    this->GlobalIDs = new INTEGER[this->N];
  }

  void Delete()
  {
    if(this->Positions != NULL)
      {
      delete [] this->Positions;
      }
    if( this->GlobalIDs != NULL )
      {
      delete [] this->GlobalIDs;
      }
    this->N = 0;
  }

  void GetBounds(REAL bounds[6])
  {
    // STEP 0: initialize bounds
    for( int i=0; i < 3; ++i )
      {
      bounds[i*2]   = std::numeric_limits<REAL>::max();
      bounds[i*2+1] = std::numeric_limits<REAL>::min();
      } // END for all dimensions

    // STEP 1: Compute bounds from particle positions
    for(INTEGER i=0; i < this->N; ++i)
      {
      for( int dim=0; dim < 3; ++dim )
        {
        if(this->Positions[i*3+dim] < bounds[dim*2])
          {
          bounds[dim*2] = this->Positions[i*3+dim];
          }
        else if(this->Positions[i*3+dim] > bounds[dim*2+1])
          {
          bounds[dim*2+1] = this->Positions[i*3+dim];
          }
        } // END for all dimensions
      } // END for all particles
  }

};

SimParticles Particles;

/**
 * @brief Computes the langrangian mesh parameters
 * @param Origin the origin of the langrangian mesh (out)
 * @param h the spacing of the langrangian mesh (out)
 * @param ext the extent of the langrangian mesh (out)
 */
void GetLangragianMeshParameters(
    REAL Origin[3], REAL h[3], INTEGER ext[6])
{
  REAL bounds[6];
  Particles.GetBounds( bounds );

  for( int i=0; i < 3; ++i )
    {
    Origin[i]  = bounds[i*2];
    h[i]       = (Parameters.rL / static_cast<REAL>(Parameters.LDIM) );
    ext[i*2]   = 0;
    ext[i*2+1] = Parameters.LDIM;
    }
}

/**
 * @brief Reads particles
 */
void ReadParticles(std::string file)
{
  Particles.Delete();

  vtkPCosmoReader *reader = vtkPCosmoReader::New();
  reader->SetRL( Parameters.rL );
  reader->SetOverlap( 0.0 );
  reader->SetFileName( file.c_str() );
  reader->Update( );

  vtkUnstructuredGrid *particles = reader->GetOutput();
  assert("pre: output particles are NULL!" && (particles != NULL) );
  vtkIntArray *tagArray =
    vtkIntArray::SafeDownCast(particles->GetPointData()->GetArray("tag"));
  int *tags = static_cast<int*>(tagArray->GetPointer(0));

  Particles.Allocate( particles->GetNumberOfPoints() );
  for( vtkIdType idx=0; idx < particles->GetNumberOfPoints(); ++idx )
    {
    Particles.Positions[idx*3]   = particles->GetPoint(idx)[0];
    Particles.Positions[idx*3+1] = particles->GetPoint(idx)[1];
    Particles.Positions[idx*3+2] = particles->GetPoint(idx)[2];

    Particles.GlobalIDs[idx] = tags[idx]-1;
    }

  reader->Delete();
}

/**
 * @brief Writes the unstructured grid to a file
 * @param m the unstructured mesh
 */
void WriteUnstructuredVTKGrid(vtkUnstructuredGrid *m, std::string file)
{
  assert("pre: m != NULL" && (m != NULL) );
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetInputData( m );
  writer->SetFileName( file.c_str() );
  writer->Update();
  writer->Delete();
}

/**
 * @brief Writes the euler mesh for the given time-step
 * @param p the structure formation probe data structure
 * @param i the given time-step
 */
void WriteEuler(
    cosmologytools::StructureFormationProbe *p, const int i)
{
  assert("pre: structure formation probe is NULL" && (p != NULL) );

  std::ostringstream oss;
  oss.clear(); oss.str("");
  oss << "euler_" << i << ".vtk";
  vtkUnstructuredGrid *eulerMesh = vtkUnstructuredGrid::New();

  std::vector<REAL> nodes;
  std::vector<INTEGER> tets;
  std::vector<REAL> vol;

  p->GetEulerMesh(nodes,tets,vol);

  assert("post: vol.size()==tets.size()/4" &&
           (vol.size()==tets.size()/4));

  INTEGER numNodes = static_cast<INTEGER>(nodes.size()/3);
  INTEGER numTets  = static_cast<INTEGER>(vol.size());

  std::cout << "Number of Euler tets: " << numTets;
  std::cout.flush();

  vtkIntArray *cellIds = vtkIntArray::New();
  cellIds->SetName("CELLID");
  cellIds->SetNumberOfComponents(1);
  cellIds->SetNumberOfTuples(numTets);

  vtkCellArray *meshElements = vtkCellArray::New();
  vtkPoints    *meshNodes    = vtkPoints::New();
  meshNodes->SetNumberOfPoints(numNodes);
  meshElements->Allocate(meshElements->EstimateSize(numTets,4));

  vtkDoubleArray *volumes = vtkDoubleArray::New();
  volumes->SetName( "Volume" );
  volumes->SetNumberOfComponents(1);
  volumes->SetNumberOfTuples(numTets);

  for(vtkIdType nodeIdx=0; nodeIdx < numNodes; ++nodeIdx)
    {
    meshNodes->SetPoint(nodeIdx,&nodes[nodeIdx*3]);
    } // END for all nodes

  vtkIdType pts[4];
  for(INTEGER tetIdx=0; tetIdx < numTets; ++tetIdx)
    {
    pts[0] = tets[tetIdx*4];
    pts[1] = tets[tetIdx*4+1];
    pts[2] = tets[tetIdx*4+2];
    pts[3] = tets[tetIdx*4+3];
    meshElements->InsertNextCell(4,pts);
    volumes->SetTuple1(tetIdx,vol[tetIdx]);
    } // END for all tets

  eulerMesh->SetPoints( meshNodes );
  meshNodes->Delete();
  eulerMesh->SetCells( VTK_TETRA, meshElements );
  meshElements->Delete();

  eulerMesh->GetCellData()->AddArray(cellIds);
  cellIds->Delete();
  eulerMesh->GetCellData()->AddArray( volumes );
  volumes->Delete();

  WriteUnstructuredVTKGrid( eulerMesh, oss.str() );
  eulerMesh->Delete();
}

/**
 * @brief Writes the faces on the caustic surfaces at the given time-step
 * @param p the structure formation probe data structure
 * @param i the given time-step
 */
void WriteCaustics(
    cosmologytools::StructureFormationProbe *p, const int i)
{
  assert("pre: structure formation probe is NULL" && (p != NULL) );

  std::ostringstream oss;
  oss.clear(); oss.str("");
  oss << "caustics_" << i << ".vtk";
  vtkUnstructuredGrid *caustics = vtkUnstructuredGrid::New();

  std::vector<REAL> nodes;
  std::vector<INTEGER> faces;
  p->ExtractCausticSurfaces(nodes,faces);

  vtkPoints *meshNodes       = vtkPoints::New();
  meshNodes->SetNumberOfPoints(nodes.size()/3);
  for(vtkIdType nodeIdx=0; nodeIdx < nodes.size()/3; ++nodeIdx)
    {
    meshNodes->SetPoint(nodeIdx,&nodes[nodeIdx*3]);
    }

  vtkCellArray *meshElements = vtkCellArray::New();
  meshElements->Allocate(meshElements->EstimateSize(faces.size()/3,3));

  vtkIdType pts[3];
  for( INTEGER cellIdx=0; cellIdx < faces.size()/3; ++cellIdx )
    {
    pts[0] = faces[ cellIdx*3  ];
    pts[1] = faces[ cellIdx*3+1];
    pts[2] = faces[ cellIdx*3+2];
    meshElements->InsertNextCell(3,pts);
    }

  caustics->SetPoints( meshNodes );
  meshNodes->Delete();
  caustics->SetCells(VTK_TRIANGLE,meshElements);
  meshElements->Delete();

  WriteUnstructuredVTKGrid(caustics, oss.str() );
  caustics->Delete();
}

/**
 * @brief Writes the uniform grid to a file
 */
void WriteUniformGrid(vtkUniformGrid *ug, std::string file)
{
  assert("pre: uniform grid should not be NULL" && (ug != NULL) );

  vtkXMLImageDataWriter *writer = vtkXMLImageDataWriter::New();

  std::ostringstream oss;
  oss.str("");
  oss.clear();
  oss << file << "." << writer->GetDefaultFileExtension();

  writer->SetFileName( oss.str().c_str() );
  writer->SetInputData( ug );
  writer->Write();
  writer->Delete();
}

/**
 * @brief Probes points on user-supplied grid (via parameters) and writes the
 * output for visualization.
 * @param p the structure formation probe data structure
 * @param i the given time-step
 */

#ifndef USEDAX

void WriteProbedGridData(
    cosmologytools::StructureFormationProbe *p, const int i)
{
  assert("pre: structure formation probe is NULL" && (p != NULL) );

  std::ostringstream oss;
  oss.clear(); oss.str("");
  oss << "probed_grid" << i << ".vtk";

  double origin[3];
  double spacing[3];
  int ext[6];

  REAL bounds[6];
  p->GetLagrangeTesselator()->GetBounds(bounds);

  for(int i=0; i < 3; ++i)
    {
    origin[i]  = bounds[i*2];
    REAL dx    = bounds[i*2+1]-bounds[i*2];
    spacing[i] = static_cast<double>(dx/Parameters.GDIM);
    ext[i*2]   = 0;
    ext[i*2+1] = Parameters.GDIM;
    }

  vtkUniformGrid *probeGrid = vtkUniformGrid::New();
  probeGrid->SetOrigin(origin);
  probeGrid->SetSpacing(spacing);
  probeGrid->SetExtent(ext);

  vtkIntArray *numberOfStreams = vtkIntArray::New();
  numberOfStreams->SetName("NumberOfStreams");
  numberOfStreams->SetNumberOfComponents( 1 );
  numberOfStreams->SetNumberOfTuples( probeGrid->GetNumberOfPoints() );
  int *nstreamPtr = numberOfStreams->GetPointer(0);

  vtkDoubleArray *rho = vtkDoubleArray::New();
  rho->SetName("rho");
  rho->SetNumberOfComponents( 1 );
  rho->SetNumberOfTuples( probeGrid->GetNumberOfPoints() );
  double *rhoPtr = rho->GetPointer(0);


  #pragma omp parallel for
  for(vtkIdType pntIdx=0; pntIdx < probeGrid->GetNumberOfPoints(); ++pntIdx )
    {
    double pnt[3];
    probeGrid->GetPoint(pntIdx,pnt);

    REAL rpnt[3];
    rpnt[0] = pnt[0];
    rpnt[1] = pnt[1];
    rpnt[2] = pnt[2];

    INTEGER nStream = 0;
    REAL lrho = 0.0;
    p->ProbePoint(rpnt,nStream,lrho);

    nstreamPtr[ pntIdx ] = static_cast<int>(nStream);
    rhoPtr[ pntIdx ]     = static_cast<double>(lrho);
    } // END for all points

  probeGrid->GetPointData()->AddArray( numberOfStreams );
  numberOfStreams->Delete();
  probeGrid->GetPointData()->AddArray( rho );
  rho->Delete();

  WriteUniformGrid( probeGrid, oss.str() );
  probeGrid->Delete();
}

#else

/**
 * @brief Probes points on user-supplied grid (via parameters) and writes the
 * output for visualization.
 * @param p the structure formation probe data structure
 * @param i the given time-step
 */
void WriteProbedGridData(
    cosmologytools::StructureFormationProbe *p, const int timestep)
{
  assert("pre: structure formation probe is NULL" && (p != NULL) );

  std::ostringstream oss;
  oss.clear(); oss.str("");
  oss << "probed_grid" << timestep << ".vtk";

  dax::Vector3 origin;
  dax::Vector3 spacing;
  dax::Extent3 ext;

  REAL bounds[6];
  p->GetLagrangeTesselator()->GetBounds(bounds);

  for(int i=0; i < 3; ++i)
    {
    origin[i]  = bounds[i*2];
    dax::Scalar dx    = bounds[i*2+1]-bounds[i*2];
    spacing[i] = static_cast<double>(dx/Parameters.GDIM);
    ext.Min[i]   = 0;
    ext.Max[i] = Parameters.GDIM;
    }

  //build the dax uniform grid from the lagrange tesselator
  dax::cont::UniformGrid<> grid;
  grid.SetOrigin(origin);
  grid.SetSpacing(spacing);
  grid.SetExtent(ext);

  std::cout << "number of points in grid: " << grid.GetNumberOfPoints() << std::endl;
  std::cout << "spacing is:  " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << std::endl;
  std::cout << "ext Min is:  " << ext.Min[0] << ", " << ext.Min[1] << ", " << ext.Min[2] << std::endl;
  std::cout << "ext Max is:  " << ext.Max[0] << ", " << ext.Max[1] << ", " << ext.Max[2] << std::endl;

  //build the dax unstructured grid from the eulermesh
  //these will be the tets we will search
  std::vector<dax::Scalar> nodes;
  std::vector<dax::Id> tetConn;
  std::vector<dax::Scalar> volumes;
  p->GetEulerMesh(nodes,tetConn,volumes);

  //to convert nodes to vec3 we need to some pointer conversion
  //this works since the memory layout of vector3 is equal to a flat
  //scalars layed out like a.x,a.y,a.z,b.x,b.y,b.z...
  dax::Vector3* startOfCoords = reinterpret_cast<dax::Vector3*>(&nodes[0]);
  dax::Id numCoordinates = nodes.size() / 3;

  typedef dax::cont::UnstructuredGrid< dax::CellTagTetrahedron > DaxTetGrid;
  DaxTetGrid tetGrid( dax::cont::make_ArrayHandle(tetConn),
                      dax::cont::make_ArrayHandle(startOfCoords ,numCoordinates ) );
  dax::cont::ArrayHandle<dax::Scalar> volumeHandle =
    dax::cont::make_ArrayHandle(volumes);

  std::cout << "EulerMesh size: " << tetGrid.GetNumberOfCells() << std::endl;


  dax::exec::TetStructure tetStructure(tetGrid,volumeHandle);

  //setup handles to the output data
  dax::cont::ArrayHandle<dax::Id> streamsHandle;
  dax::cont::ArrayHandle<dax::Scalar> rhoHandle;

  dax::cont::Scheduler<> scheduler;

  dax::cont::Timer<> timer;
  scheduler.Invoke(dax::worklet::MapTetsToPoints(),
                   grid.GetPointCoordinates(),
                   tetStructure,
                   streamsHandle,
                   rhoHandle);

 std::cout << "Time to process one point: " << timer.GetElapsedTime() << std::endl;

//todo implement writing from dax
// WriteUniformGrid( grid, streamsHandle, rhoHandle, oss.str() );
}

#endif

//=============================================================================

/**
 * @brief Program main
 * @param argc argument counter
 * @param argv argument vector
 * @return rc return code
 */
int main(int argc, char **argv)
{
  // STEP 0: Read program arguments
  if( argc != 6)
    {
    std::cerr << "USAGE: ./probe <basedir> <rL> <ndim> <ndim2> <fringe>\n";
    }

  Parameters.BaseDir = std::string(argv[1]);
  Parameters.rL      = atof(argv[2]);
  Parameters.LDIM    = atoi(argv[3]);
  Parameters.GDIM    = atoi(argv[4]);
  Parameters.Fringe  = atoi(argv[5]);
  Parameters.Println();

  // STEP 1: Setup structure formation probe
  cosmologytools::StructureFormationProbe *Probe =
      new cosmologytools::StructureFormationProbe();
  Probe->SetFringe( Parameters.Fringe );

  // STEP 2: For each timestep
  for( int i=0; i < NSTEPS; ++i )
    {
    std::string file = Parameters.BaseDir + files[i];
    std::cout << "==================================================\n";
    std::cout << "tstep=" << i << " reading file " << file << "...";
    ReadParticles( file );
    std::cout << "[DONE]\n";
    std::cout.flush();

    std::cout << "Total Number of particles: " << Particles.N << std::endl;
    std::cout.flush();
    if( Particles.N == 0 )
      {
      continue;
      }

    // STEP 3: Build langrangian mesh if this is the first time
    if( i == 0 )
      {
      REAL Origin[3]; REAL h[3]; INTEGER ext[6];
      GetLangragianMeshParameters(Origin,h,ext);

      std::cout << "Building Langrangian mesh...\n";
      std::cout << "ORIGIN=" << Origin[0] << " ";
      std::cout << Origin[1] << " " << Origin[2] << std::endl;
      std::cout << "h=" << h[0] << " " << h[1] << " " << h[2] << std::endl;
      std::cout << "EXTENT=[ ";
      for( int dim=0; dim < 3; ++dim )
        {
        std::cout << ext[dim*2] << " " << ext[dim*2+1] << " ";
        }
      std::cout << "]\n";
      std::cout.flush();

      Probe->BuildLangrangianMesh(Origin,h,ext);
      std::cout << "Number of Langrangian tets: ";
      std::cout << Probe->GetLagrangeTesselator()->GetNumTets() << std::endl;
      std::cout.flush();
      }

    // STEP 4: Set particles for this time-step
    std::cout << "Set Particles...";
    std::cout.flush();
    Probe->SetParticles(
        Particles.Positions,Particles.GlobalIDs,Particles.N);
    std::cout << "[DONE]\n";
    std::cout.flush();

    // STEP 5: Build euler mesh
    std::cout << "Build euler mesh...";
    std::cout.flush();
    Probe->BuildEulerMesh();
    std::cout << "[DONE]\n";
    std::cout.flush();

    // STEP 6: Write euler mesh
    std::cout << "Write Euler mesh...";
    std::cout.flush();
    WriteEuler(Probe, i);
    std::cout << "[DONE]\n";
    std::cout.flush();

    // STEP 7: Write caustics
    std::cout << "Extract caustic surfaces...";
    std::cout.flush();
    WriteCaustics(Probe,i);
    std::cout << "[DONE]\n";
    std::cout.flush();

    // STEP 8: Write probed data
    std::cout << "Write Probed grid...";
    std::cout.flush();
    WriteProbedGridData(Probe,i);
    std::cout << "[DONE]\n";
    std::cout.flush();

    // STEP 9: Statistics:
    std::cout << "NumPoints probed=" << Probe->GetNumPointsProbed();
    std::cout << std::endl;
    std::cout << "NumTets checked=" << Probe->GetNumTetsChecked();
    std::cout << std::endl;
    std::cout << "Avg. NumTets per point=";
    std::cout << Probe->GetNumTetsChecked()/Probe->GetNumPointsProbed();
    std::cout << std::endl;
    }

  delete Probe;
  return 0;
}

