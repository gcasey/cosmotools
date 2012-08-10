/*=========================================================================

 Program:   Visualization Toolkit
 Module:    vtkPStructureFormationProbe.h

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/
// .NAME vtkPStructureFormationProbe.h -- Structure formation filter
//
// .SECTION Description
//  A concrete instance of vtkMultiBlockDataSetAlgorithm that implements a
//  tesselation-based approach for probing the structure of the universe. Namely,
//  the filter implements functionality for:
//  <ol>
//    <li>Finding streams & extracting associated regions, e.g., filaments</li>
//    <li>Extracting caustic surfaces</li>
//  </ol>
//  The input of the filter is a set of particles and the output consists of the
//  resulting tesselation in langrangian (initial) or euler space.
#ifndef VTKPSTRUCTUREFORMATIONPROBE_H_
#define VTKPSTRUCTUREFORMATIONPROBE_H_

#include "vtkUnstructuredGridAlgorithm.h"

// sfprobe includes
#include "CosmologyToolsMacros.h"
#include "StructureFormationProbe.h"

// C/C++ includes
#include <map>

// Forward declarations
class vtkInformation;
class vtkInformationVector;
class vtkMultiProcessController;
class vtkUniformGrid;
class vtkUnstructuredGrid;

class VTK_EXPORT vtkPStructureFormationProbe :
  public vtkUnstructuredGridAlgorithm
{
public:

  enum {
    LANGRANGE = 0,
    EULER     = 1
  } Space;

  static vtkPStructureFormationProbe *New();
  vtkTypeMacro(vtkPStructureFormationProbe,vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set the extent imin,imax,jmin,jmax,kmin,kmax
  vtkSetVector6Macro(Extent,int);
  vtkGetVector6Macro(Extent,int);

  // Description:
  // Set the domain space, either LANGRANGE or EULER.
  vtkSetMacro(DomainSpace,int);
  vtkGetMacro(DomainSpace,int);

  // Description:
  // Sets periodic boundary fringe. Periodic boundary fringe
  vtkSetMacro(Fringe,int);
  vtkGetMacro(Fringe,int);

  // Description:
  // Set/Get macro for the controller.
  vtkSetMacro(Controller,vtkMultiProcessController*);
  vtkGetMacro(Controller,vtkMultiProcessController*);

  // Description:
  // Set/Get whether or not to shift the global numbering to start from 0.
  vtkSetMacro(ShiftGlobalNumberingToZero,int);
  vtkGetMacro(ShiftGlobalNumberingToZero,int);

  // Description:
  // Set/Get whether to probe the langrangian mesh.
  vtkSetMacro(ProbeLangrangianMesh,int);
  vtkGetMacro(ProbeLangrangianMesh,int);

protected:
  vtkPStructureFormationProbe();
  virtual ~vtkPStructureFormationProbe();

  // Description:
  // Standard pipeline operations
  virtual int RequestData(
      vtkInformation *request,vtkInformationVector **inputVector,
      vtkInformationVector *outputVector );
  virtual int FillInputPortInformation(int port, vtkInformation *info);
  virtual int FillOutputPortInformation( int port, vtkInformation *info );

  // Description:
  // Extracts the input particle dataset from the unstructured grid.
  void ExtractInputDeck( vtkUnstructuredGrid *particles );

  // Description:
  // Builds a langrangian tesselation
  void BuildLangrangeTesselation();

  // Description:
  // Constructs the output mesh. Depending on what DomainSpace the user has
  // request, this method will construct the output mesh in lagrangian or
  // euler space.
  void ConstructOutputMesh(vtkUnstructuredGrid *tesselation);

  // Description:
  // Gets the euler mesh
  void GetEulerMesh(vtkUnstructuredGrid *eulerMesh);

  // Description:
  // Gets the langrangian mesh
  void GetLangrangianMesh(vtkUnstructuredGrid *langrangeMesh);

  // Description:
  // Inspects the euler-grid and extracts the faces on caustic surfaces.
  // Caustic surfaces are defined as by the (triangular) faces of two adjacent
  // tetrahedra that have a volume with opposite signs. The caustic surfaces
  // can be identified only in euler space, however, the output is given either
  // in euler or langrangian space depending on the user-supplied space
  // parameter.
  void ExtractCausticSurfaces(
        vtkUnstructuredGrid *particles,
        vtkUnstructuredGrid *tesselation,
        vtkUnstructuredGrid *caustics);

  // Description:
  // Probes the langrangian grid on the euler mesh and computes the number of
  // streams and local density at each grid point.
  void ProbeGrid(vtkUniformGrid *probeGrid);

  // User-supplied parameters
  int DomainSpace;
  int ShiftGlobalNumberingToZero;
  int Extent[6];
  int Fringe;
  int ProbeLangrangianMesh;

  // Computed ivars
  double Bounds[6];
  REAL *Particles;
  INTEGER *GlobalIds;
  INTEGER N;

  vtkMultiProcessController *Controller;
  cosmologytools::StructureFormationProbe *SFProbe;
private:
  vtkPStructureFormationProbe(const vtkPStructureFormationProbe&); // Not implemented
  void operator=(const vtkPStructureFormationProbe&); // Not implemented
};

#endif /* VTKPSTRUCTUREFORMATIONPROBE_H_ */
