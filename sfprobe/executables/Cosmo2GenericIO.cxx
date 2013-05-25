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


/**
 * @brief Program main
 * @param argc argument counter
 * @param argv argument vector
 * @return rc return code
 */
int main(int argc, char **argv)
{
  return 0;
}

/**
 * @brief Reads particles
 */
void ReadParticles(std::string file, vtkUnstructuredGrid *partciles)
{
  Particles.Delete();

  vtkPCosmoReader *reader = vtkPCosmoReader::New();
  reader->SetRL( Parameters.rL );
  reader->SetOverlap( 0.0 );
  reader->SetFileName( file.c_str() );
  reader->Update( );

  particles = reader->GetOutput();
  reader->Delete();
}
