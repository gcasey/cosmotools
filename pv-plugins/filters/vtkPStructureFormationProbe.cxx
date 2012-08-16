#include "vtkPStructureFormationProbe.h"

// VTK includes
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDoubleArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkMultiProcessController.h"
#include "vtkObjectFactory.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkStructuredData.h"
#include "vtkUniformGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkMultiBlockDataSet.h"


// sfprobe includes
#include "ExtentUtilities.h"
#include "TetrahedronUtilities.h"
#include "LangrangianTesselator.h"

// C/C++ includes
#include <cassert>

vtkStandardNewMacro(vtkPStructureFormationProbe);

vtkPStructureFormationProbe::vtkPStructureFormationProbe()
{
   for( int i=0; i < 6; ++i )
    {
    this->Extent[i] = 0;
    }

  this->DomainSpace                = EULER;
  this->SFProbe                    = NULL;
  this->Particles                  = NULL;
  this->GlobalIds                  = NULL;
  this->N                          = 0;
  this->ShiftGlobalNumberingToZero = 1;
  this->Fringe                     = 1;
  this->CurrentFringe              = -1;
  this->ProbeGrid                  = 1;
  this->CurrentTimeStep            = -1;

  this->EulerGrid      = NULL;
  this->CausticSurface = NULL;
  this->AuxiliaryGrid  = NULL;

  this->Controller = vtkMultiProcessController::GetGlobalController();
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(3);
}

//------------------------------------------------------------------------------
vtkPStructureFormationProbe::~vtkPStructureFormationProbe()
{
  if( this->SFProbe != NULL )
    {
    delete this->SFProbe;
    }

  if( this->Particles != NULL )
    {
    delete [] this->Particles;
    }

  if( this->GlobalIds != NULL )
    {
    delete [] this->GlobalIds;
    }

  if( this->EulerGrid != NULL )
    {
    this->EulerGrid->Delete();
    }

  if( this->CausticSurface != NULL )
    {
    this->CausticSurface->Delete();
    }

  if( this->AuxiliaryGrid != NULL )
    {
    this->AuxiliaryGrid->Delete();
    }
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::PrintSelf(
          std::ostream &os, vtkIndent indent)
{
  this->Superclass::PrintSelf( os, indent );
}

//------------------------------------------------------------------------------
int vtkPStructureFormationProbe::FillInputPortInformation(
                                  int vtkNotUsed(port), vtkInformation *info)
{
  assert( "pre: info is NULL" && (info != NULL) );
  info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(),"vtkUnstructuredGrid");
  return 1;
}

//------------------------------------------------------------------------------
int vtkPStructureFormationProbe::FillOutputPortInformation(
                                  int port, vtkInformation *info)
{
  assert( "pre: info is NULL" && (info != NULL) );
  switch( port )
    {
    case 0:
    case 1:
      info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkUnstructuredGrid");
      break;
    case 2:
      info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkMultiBlockDataSet");
      break;
    default:
      vtkErrorMacro("output port index out-of-bounds!");
    }
  return 1;
}

//------------------------------------------------------------------------------
int vtkPStructureFormationProbe::RequestData(
    vtkInformation *request,vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  std::cout << "===\n";
  std::cout.flush();

  // STEP 0: Get input object
  vtkInformation *input = inputVector[0]->GetInformationObject( 0 );
  assert("pre: input information object is NULL" && (input != NULL) );
  vtkUnstructuredGrid *particles =
      vtkUnstructuredGrid::SafeDownCast(
          input->Get( vtkDataObject::DATA_OBJECT() ) );
  assert("pre: input particles is NULL!" && (particles != NULL) );

  // STEP 1: Get current time-step
  int timeStep  = 0;
  if(input->Has(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()))
    {
    timeStep = static_cast<int>(
      input->Get(vtkStreamingDemandDrivenPipeline::UPDATE_TIME_STEP()));
    }
  std::cout << "TimeStep=" << timeStep << std::endl;
  std::cout.flush();

  // STEP 2: Get 1st output object (tesselation)
  vtkInformation *output1 = outputVector->GetInformationObject( 0 );
  assert( "pre: output information object is NULL" && (output1 != NULL) );
  vtkUnstructuredGrid *tesselation =
      vtkUnstructuredGrid::SafeDownCast(
          output1->Get( vtkDataObject::DATA_OBJECT() ) );

  // STEP 3: Get 2nd output object (surface caustics)
  vtkInformation *output2 = outputVector->GetInformationObject( 1 );
  assert( "pre: output information object is NULL" && (output2 != NULL) );
  vtkUnstructuredGrid *caustics =
      vtkUnstructuredGrid::SafeDownCast(
          output2->Get( vtkDataObject::DATA_OBJECT() ) );

  // STEP 4: Get 3rd output (probe grid)
  vtkInformation *output3 = outputVector->GetInformationObject( 2 );
  assert( "pre: output information object is NULL" && (output3 != NULL) );
  vtkMultiBlockDataSet *probedGrid =
      vtkMultiBlockDataSet::SafeDownCast(
            output3->Get(vtkDataObject::DATA_OBJECT() ) );
  // NOTE: The number of blocks should be the same as the number of ranks when
  // we eventually run in parallel
  probedGrid->SetNumberOfBlocks( 1 );

  // Flag used to indicate that the remaining steps of the algorithm must
  // re-execute
  bool reExecute = false;

  // STEP 5: Check if we need to re-execute for the given time-step
  if( this->NewTimeStep(timeStep) )
    {
    this->CurrentTimeStep = timeStep;
    reExecute = true;
    }

  // STEP 6: Extract input deck to the structure formation probe code
  this->ExtractInputDeck( particles );
  if( this->SFProbe == NULL )
    {
    this->SFProbe = new cosmologytools::StructureFormationProbe();
    }

  // STEP 7: Build Langrange tesselation
  // This should execute iff it is the first time-step, or if the langrangian
  // grid parameters have changed
  if( this->SFProbe->GetLangrangeTesselator() == NULL ||
       this->ExtentsAreDifferent(
         this->Extent,this->CurrentLangrangianExtent) )
    {
    // Update current langrangian extent
    this->SetCurrentLangrangianExtent( this->Extent );

    std::cout << "Building Langrangian mesh...";
    this->BuildLangrangeTesselation();
    std::cout << "[DONE]\n";
    std::cout.flush();

    // Force down-stream calculations to execute
    reExecute = true;
    }

  // STEP 8: Construct the output mesh
  // This should execute iff the fringe has changed, or something changed
  // upstream, i.e., the langrangian mesh
  if( reExecute || this->FringeHasChanged() )
    {
    if( this->EulerGrid != NULL )
      {
      this->EulerGrid->Delete();
      }
    this->EulerGrid = vtkUnstructuredGrid::New();

    std::cout << "Build output tesselation mesh...";
    this->CurrentFringe = this->Fringe;
    this->SFProbe->SetFringe( this->Fringe );
    this->ConstructOutputMesh(this->EulerGrid);
    std::cout << "[DONE]\n";
    std::cout.flush();
    // Force down-stream calculations to execute
    reExecute = true;
    }

  // STEP 9: Extract the caustic surfaces
  if( reExecute )
    {
    // Note, it only makese sent to extract caustic faces in euler space
    if( this->DomainSpace == EULER )
      {
      if( this->CausticSurface != NULL )
        {
        this->CausticSurface->Delete();
        }
      this->CausticSurface = vtkUnstructuredGrid::New();

      std::cout << "Extracting caustic surfaces...";
      this->ExtractCausticSurfaces(particles,
          this->EulerGrid,this->CausticSurface);
      std::cout << "[DONE]\n";
      std::cout.flush();
      } // END if in Euler space
    }


  // STEP 10: Probe grid
  if( reExecute || this->ExtentsAreDifferent(
                      this->ProbeGridExtent,this->CurrentProbeGridExtent))
    {

    if( this->AuxiliaryGrid != NULL )
      {
      this->AuxiliaryGrid->Delete();
      }
    this->AuxiliaryGrid = vtkUniformGrid::New();

    // Update current probe grid extent
    this->SetCurrentProbeGridExtent( this->ProbeGridExtent );

    // Note, it only makes sense to probe the grid in euler space
    if( this->DomainSpace == EULER )
      {
      std::cout << "Probing uniform grid...";
      this->ProbeUniformGrid( this->AuxiliaryGrid );
      std::cout << "[DONE]\n";
      std::cout.flush();
      } // END if in Euler space
    }

  // Shallow-copy output objects
  tesselation->ShallowCopy( this->EulerGrid );
  caustics->ShallowCopy( this->CausticSurface );

  // Set probe grid output
  probedGrid->SetBlock(0,this->AuxiliaryGrid);

  this->Controller->Barrier();
  return 1;
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::ExtractInputDeck(
        vtkUnstructuredGrid *particles )
{
  assert("pre: input particles dataset is NULL" && (particles != NULL) );
  assert("pre: no 'tag' array associated with the input particle dataset!"
          && (particles->GetPointData()->HasArray("tag")) );

  vtkIntArray *tagArray =
      vtkIntArray::SafeDownCast(particles->GetPointData()->GetArray("tag"));
  int *tags = static_cast<int*>(tagArray->GetPointer(0));

  this->N = particles->GetNumberOfPoints();
  if( this->Particles != NULL )
    {
    delete [] this->Particles;
    }
  this->Particles = new REAL[3*this->N];

  if( this->GlobalIds != NULL )
    {
    delete [] this->GlobalIds;
    }
  this->GlobalIds = new INTEGER[this->N];

  particles->GetBounds( this->Bounds );

  for(vtkIdType idx=0; idx < particles->GetNumberOfPoints(); ++idx)
    {
    this->Particles[idx*3]   = particles->GetPoint(idx)[0];
    this->Particles[idx*3+1] = particles->GetPoint(idx)[1];
    this->Particles[idx*3+2] = particles->GetPoint(idx)[2];

    this->GlobalIds[idx] = tags[idx];

    if( this->ShiftGlobalNumberingToZero==1 )
      {
      this->GlobalIds[idx]--;
      }
    } // END for all particles
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::GetLangrangianMesh(
    vtkUnstructuredGrid *langrangeMesh)
{
  assert("pre: Langrangian tesselator has not been allocated" &&
         (this->SFProbe->GetLangrangeTesselator() != NULL) );
  assert("pre: input langrange mesh is NULL" && (langrangeMesh != NULL) );
  assert("pre: Number of langrange tets > 0" &&
         (this->SFProbe->GetLangrangeTesselator()->GetNumTets() > 0) );

  vtkIntArray *cellIds = vtkIntArray::New();
  cellIds->SetName("CELLID");
  cellIds->SetNumberOfComponents(1);
  cellIds->SetNumberOfTuples(
      this->SFProbe->GetLangrangeTesselator()->GetNumTets());

  vtkCellArray *meshElements = vtkCellArray::New();
  vtkPoints    *meshNodes    = vtkPoints::New();
  meshNodes->SetNumberOfPoints(
      vtkStructuredData::GetNumberOfNodes(this->Extent));
  meshElements->Allocate(
      meshElements->EstimateSize(
          this->SFProbe->GetLangrangeTesselator()->GetNumTets(),4));

  vtkDoubleArray *volumes = NULL;
  volumes = vtkDoubleArray::New();
  volumes->SetName( "Volume" );
  volumes->SetNumberOfComponents(1);
  volumes->SetNumberOfTuples(
      this->SFProbe->GetLangrangeTesselator()->GetNumTets());


  vtkIdList *tetElement = vtkIdList::New();
  tetElement->SetNumberOfIds(4);
  for(INTEGER idx=0;
      idx < this->SFProbe->GetLangrangeTesselator()->GetNumTets(); ++idx )
    {
    INTEGER tet[4];
    REAL v0[3]; REAL v1[3]; REAL v2[3]; REAL v3[3];

    this->SFProbe->GetLangrangeTesselator()->
                      GetLangrangianTet(idx,v0,v1,v2,v3,tet);
    assert("pre: volumes arrays has not been allocated" &&
           (volumes != NULL));
    REAL vol=
        cosmologytools::TetrahedronUtilities::ComputeVolume(v0,v1,v2,v3);
    assert("pre: tet in langrangian space with negative volume!" &&
           (vol > 0.0) );
    volumes->SetTuple1(idx,static_cast<double>(vol));


    meshNodes->SetPoint(tet[0],v0);
    meshNodes->SetPoint(tet[1],v1);
    meshNodes->SetPoint(tet[2],v2);
    meshNodes->SetPoint(tet[3],v3);

    for( int i=0; i < 4; ++i )
      {
      tetElement->SetId(i,tet[i]);
      }

    cellIds->SetValue(idx,idx);
    meshElements->InsertNextCell(tetElement);
    } // END for all tets
  tetElement->Delete();

  langrangeMesh->SetPoints( meshNodes );
  meshNodes->Delete();
  langrangeMesh->SetCells( VTK_TETRA, meshElements );
  meshElements->Delete();

  langrangeMesh->GetCellData()->AddArray(cellIds);
  cellIds->Delete();
  langrangeMesh->GetCellData()->AddArray( volumes );
  volumes->Delete();
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::ProbeUniformGrid(
        vtkUniformGrid *probeGrid)
{
  assert("pre: Structure formation probe is NULL" && (this->SFProbe != NULL));
  assert("pre: probe grid should not be NULL" && (probeGrid != NULL) );

  if( this->ProbeGrid != 1 )
    {
    return;
    }

  // STEP 0: Get langrangian grid parameters, origin, spacing etc.
  double origin[3];
  double spacing[3];
  int ext[6];

  cosmologytools::LangrangianTesselator *langrangianMesh =
                            this->SFProbe->GetLangrangeTesselator();
  assert("pre: langrangian mesh should not be NULL!" &&
          (langrangianMesh != NULL) );

  REAL bounds[6];
  langrangianMesh->GetBounds(bounds);

  for( int i=0; i < 3; ++i )
    {
    origin[i] = bounds[i*2];
    REAL dx   = bounds[i*2+1]-bounds[i*2];
    ext[ i*2   ] = this->ProbeGridExtent[i*2];
    ext[ i*2+1 ] = this->ProbeGridExtent[i*2+1];
    int ndim     = ext[i*2+1]-ext[i*2]+1;
    spacing[i]   = static_cast<double>(dx/static_cast<double>(ndim));
    }

  // STEP 1: Construct probe grid
  probeGrid->SetOrigin(origin);
  probeGrid->SetSpacing(spacing);
  probeGrid->SetExtent(ext);

  // STEP 2: Create arrays to store the number of streams and the local
  // density rho
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

  // STEP 3: Loop through each grid point and probe it on the euler mesh
  vtkIdType pntIdx = 0;
  for( ; pntIdx < probeGrid->GetNumberOfPoints(); ++pntIdx )
    {
    double pnt[3];
    probeGrid->GetPoint(pntIdx,pnt);

    REAL rpnt[3];
    rpnt[0] = pnt[0];
    rpnt[1] = pnt[1];
    rpnt[2] = pnt[2];

    INTEGER nStream = 0;
    REAL lrho = 0.0;
    this->SFProbe->ProbePoint(rpnt,nStream,lrho);

    nstreamPtr[ pntIdx ] = static_cast<int>(nStream);
    rhoPtr[ pntIdx ]     = static_cast<double>(lrho);
    } // END for all points

  // STEP 4: Add probed data arrays to grid
  probeGrid->GetPointData()->AddArray( numberOfStreams );
  numberOfStreams->Delete();
  probeGrid->GetPointData()->AddArray( rho );
  rho->Delete();
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::ConstructOutputMesh(
    vtkUnstructuredGrid *tess)
{
  assert("pre: input tesselation object is NULL" && (tess != NULL) );

  switch( this->DomainSpace )
    {
    case EULER:
      this->GetEulerMesh(tess);
      break;
    case LANGRANGE:
      this->GetLangrangianMesh(tess);
      break;
    default:
      vtkErrorMacro("Undefined domain space");
    }
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::GetEulerMesh(
          vtkUnstructuredGrid *eulerMesh)
{
  this->SFProbe->SetParticles(this->Particles,this->GlobalIds,this->N);
  this->SFProbe->BuildEulerMesh();


  std::vector<REAL> nodes;
  std::vector<INTEGER> tets;
  std::vector<REAL> vol;

  this->SFProbe->GetEulerMesh(nodes,tets,vol);

  assert("post: vol.size()==tets.size()/4" &&
           (vol.size()==tets.size()/4));

  INTEGER numNodes = static_cast<INTEGER>(nodes.size()/3);
  INTEGER numTets  = static_cast<INTEGER>(vol.size());

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
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::ExtractCausticSurfaces(
        vtkUnstructuredGrid *particles,
        vtkUnstructuredGrid *tess,
        vtkUnstructuredGrid *caustics)
{
  assert("pre: input particles is NULL" && (particles != NULL) );
  assert("pre: input tesselation object is NULL" && (tess != NULL) );
  assert("pre: input caustics object is NULL" && (caustics != NULL) );

  std::vector<REAL> nodes;
  std::vector<INTEGER> faces;
  this->SFProbe->ExtractCausticSurfaces(nodes,faces);

  vtkPoints *meshNodes       = vtkPoints::New();
  meshNodes->SetNumberOfPoints(nodes.size());
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
}

//------------------------------------------------------------------------------
void vtkPStructureFormationProbe::BuildLangrangeTesselation()
{
  REAL origin[3];
  REAL h[3];
  INTEGER ext[6];

  int dims[3];
  vtkStructuredData::GetDimensionsFromExtent(this->Extent,dims);
  //Convert to internal types that the CosmologyTools library was compiled with
  for( int i=0; i < 3; ++i )
    {
    origin[ i ] = this->Bounds[i*2];
    h[ i ]      = (this->Bounds[i*2+1]-this->Bounds[i*2]) /
                   static_cast<double>(dims[i]);
    ext[i*2]    = static_cast<INTEGER>( this->Extent[i*2] );
    ext[i*2+1]  = static_cast<INTEGER>( this->Extent[i*2+1]);
    }

  this->SFProbe->BuildLangrangianMesh(origin,h,ext);
}
