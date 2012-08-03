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
#include "vtkStructuredData.h"
#include "vtkUnstructuredGrid.h"


// sfprobe includes
#include "ExtentUtilities.h"
#include "TetrahedronUtilities.h"

// C/C++ includes
#include <cassert>

vtkStandardNewMacro(vtkPStructureFormationProbe);

vtkPStructureFormationProbe::vtkPStructureFormationProbe()
{
  this->DomainSpace = LANGRANGE;
  for( int i=0; i < 6; ++i )
    {
    this->Extent[i] = 0;
    }

  this->SFProbe                    = NULL;
  this->Particles                  = NULL;
  this->GlobalIds                  = NULL;
  this->N                          = 0;
  this->ShiftGlobalNumberingToZero = 1;
  this->Fringe                     = 1;

  this->Controller = vtkMultiProcessController::GetGlobalController();
  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(2);
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
    delete this->Particles;
    }

  if( this->GlobalIds != NULL )
    {
    delete this->GlobalIds;
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
                                  int vtkNotUsed(port), vtkInformation *info)
{
  assert( "pre: info is NULL" && (info != NULL) );
  info->Set(vtkDataObject::DATA_TYPE_NAME(),"vtkUnstructuredGrid");
  return 1;
}

//------------------------------------------------------------------------------
int vtkPStructureFormationProbe::RequestData(
    vtkInformation *request,vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{
  // STEP 0: Get input object
  vtkInformation *input = inputVector[0]->GetInformationObject( 0 );
  assert("pre: input information object is NULL" && (input != NULL) );
  vtkUnstructuredGrid *particles =
      vtkUnstructuredGrid::SafeDownCast(
          input->Get( vtkDataObject::DATA_OBJECT() ) );
  assert("pre: input particles is NULL!" && (particles != NULL) );

  // STEP 1: Get 1st output object (tesselation)
  vtkInformation *output1 = outputVector->GetInformationObject( 0 );
  assert( "pre: output information object is NULL" && (output1 != NULL) );
  vtkUnstructuredGrid *tesselation =
      vtkUnstructuredGrid::SafeDownCast(
          output1->Get( vtkDataObject::DATA_OBJECT() ) );

  // STEP 2: Get 2nd output object (surface caustics)
  vtkInformation *output2 = outputVector->GetInformationObject( 1 );
  assert( "pre: output information object is NULL" && (output2 != NULL) );
  vtkUnstructuredGrid *caustics =
      vtkUnstructuredGrid::SafeDownCast(
          output2->Get( vtkDataObject::DATA_OBJECT() ) );

  // STEP 3: Extract input deck to the structure formation probe code
  this->ExtractInputDeck( particles );
  if( this->SFProbe == NULL )
    {
    this->SFProbe = new cosmologytools::StructureFormationProbe();
    this->SFProbe->SetFringe( this->Fringe );
    }

  // STEP 4: Build Langrange tesselation
  this->BuildLangrangeTesselation();

  // STEP 5: Construct the output mesh
  this->ConstructOutputMesh(tesselation);

  // STEP 6: Extract the caustic surfaces
  this->ExtractCausticSurfaces(particles, tesselation, caustics);

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
  std::cout << "Extracting caustic surfaces....";
  this->SFProbe->ExtractCausticSurfaces(nodes,faces);
  std::cout << "[DONE]\n";
  std::cout.flush();

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
