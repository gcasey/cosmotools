/*=========================================================================

  Program:   Visualization Toolkit
  Module:    mtreeviz.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

// C/C++ includes
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cassert>

// STL includes
#include <vector>
#include <map>

// VTK includes
#include "vtkActor.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkGraphLayout.h"
#include "vtkGraphLayoutView.h"
#include "vtkGraphToPolyData.h"
#include "vtkIdList.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkPointData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataWriter.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTree.h"
#include "vtkTreeLayoutStrategy.h"
#include "vtkUnstructuredGrid.h"
#include "vtkUnstructuredGridWriter.h"

struct TreeNodeInfo {

  // Data read from file
  float HaloNumber;
  float OriginalHaloIdx;
  float DescendantNumber;
  float TimeStep;
  float Volume;
  float PosX;
  float PosY;
  float PosZ;
  float VelX;
  float VelY;
  float VelZ;

  /**
   * @brief Returns the position for the graph node
   * @param pos the graph node position (out)
   */
  void GetGraphNodePosition(double pos[3])
    {
    pos[0] = static_cast<double>(this->OriginalHaloIdx);
    pos[1] = static_cast<double>(this->TimeStep);
    pos[2] = 0.0;
    }

  /**
   * @return HashCode for
   */
 int GetHashCode()
   {
   return( static_cast<int>(this->HaloNumber) );
   }

 /**
  * @brief Prints the info to the console.
  * @note For debugging.
  */
 void Println()
   {
   std::cout << "=======\n";
   std::cout << "HALO ID: "         << this->HaloNumber       << std::endl;
   std::cout << "ORIGINAL HALO ID:" << this->OriginalHaloIdx << std::endl;
   std::cout << "DESCENDANT: "      << this->DescendantNumber << std::endl;
   std::cout << "Time-step: "       << this->TimeStep         << std::endl;
   std::cout << "Volume: "          << this->Volume           << std::endl;
   std::cout << "PosX: "            << this->PosX             << std::endl;
   std::cout << "PosY: "            << this->PosY             << std::endl;
   std::cout << "PosZ: "            << this->PosZ             << std::endl;
   std::cout << "VelX: "            << this->VelX             << std::endl;
   std::cout << "VelY: "            << this->VelY             << std::endl;
   std::cout << "VelZ: "            << this->VelZ             << std::endl;
   std::cout << "hashcode: "        << this->GetHashCode()    << std::endl;
   std::cout.flush();
   }

};

// Stores all the nodes in the graph
std::vector< TreeNodeInfo > GraphNodes;

vtkUnstructuredGrid *graphGrid;


/**
 * @brief Reads data into
 * @param fileName
 * @return
 */
bool ReadData( std::string fileName )
{
  std::ifstream ifs(fileName.c_str(),std::ios::in|std::ios::binary);
  assert("pre: cannot open input file!" && ifs.is_open() );

  if( !ifs.is_open() )
    {
    return false;
    }

  TreeNodeInfo tnode;
  while( !ifs.eof() )
    {
    ifs.read((char*)(&tnode.HaloNumber),sizeof(float));
    ifs.read((char*)(&tnode.OriginalHaloIdx),sizeof(float));
    ifs.read((char*)(&tnode.DescendantNumber),sizeof(float));
    ifs.read((char*)(&tnode.TimeStep),sizeof(float));
    ifs.read((char*)(&tnode.Volume),sizeof(float));
    ifs.read((char*)(&tnode.PosX),sizeof(float));
    ifs.read((char*)(&tnode.PosY),sizeof(float));
    ifs.read((char*)(&tnode.PosZ),sizeof(float));
    ifs.read((char*)(&tnode.VelX),sizeof(float));
    ifs.read((char*)(&tnode.VelY),sizeof(float));
    ifs.read((char*)(&tnode.VelZ),sizeof(float));

    if( ifs.eof() )
      {
      break;
      }

    tnode.Println();
    GraphNodes.push_back( tnode );

    } // END while !ifs.eof()

  std::cout << "Number of nodes: " << GraphNodes.size() << std::endl;
  std::cout.flush();

  ifs.close();
  return true;
}

/**
 * @brief Builds the graph.
 * @param g the graph data-structure to be populated.
 */
void BuildGraph(vtkMutableDirectedGraph *g)
{
  assert("pre: graph should not be NULL!" && (g != NULL) );

  graphGrid = vtkUnstructuredGrid::New();
  vtkCellArray *meshLines = vtkCellArray::New();
  meshLines->Allocate(meshLines->EstimateSize(GraphNodes.size()-1,2));


  vtkDoubleArray *volumes = vtkDoubleArray::New();
  volumes->SetName("HaloDensity");
  volumes->SetNumberOfComponents( 1 );
  volumes->SetNumberOfTuples( GraphNodes.size() );

  vtkDoubleArray *velocities = vtkDoubleArray::New();
  velocities->SetName("Velocities");
  velocities->SetNumberOfComponents( 3 );
  velocities->SetNumberOfTuples( GraphNodes.size() );

  vtkPoints *points = vtkPoints::New();

  double *velptr = static_cast<double*>(velocities->GetPointer(0));
  double *vptr   = static_cast<double*>(volumes->GetPointer(0));

  // STEP 0: Insert nodes
  double pos[3];
  for( int i=0; i < GraphNodes.size(); ++i )
    {
    GraphNodes[i].GetGraphNodePosition(pos);
    vptr[ i ] = GraphNodes[ i ].Volume;

    velptr[ i*3   ] = GraphNodes[ i ].VelX;
    velptr[ i*3+1 ] = GraphNodes[ i ].VelY;
    velptr[ i*3+2 ] = GraphNodes[ i ].VelZ;

    vtkIdType idx = g->AddVertex();
    points->InsertNextPoint( pos );
    } // END for all graph nodes
  g->SetPoints( points );
  graphGrid->SetPoints( points );
  points->Delete();

  g->GetVertexData()->AddArray( volumes );
  graphGrid->GetPointData()->AddArray( volumes );
  volumes->Delete();

  g->GetVertexData()->AddArray( velocities );
  graphGrid->GetPointData()->AddArray( velocities );
  velocities->Delete();


  // STEP 1: Create edges
  vtkIdList *edgeNodes = vtkIdList::New();
  edgeNodes->SetNumberOfIds(2);
  for( int i=0; i < GraphNodes.size(); ++i )
    {
    if( static_cast<int>(GraphNodes[i].DescendantNumber) != -1)
      {
      g->AddEdge( i, GraphNodes[i].DescendantNumber);
      edgeNodes->SetId(0,i);
      edgeNodes->SetId(1,GraphNodes[i].DescendantNumber);
      meshLines->InsertNextCell( edgeNodes );
      }
    } // END for all graph nodes
  edgeNodes->Delete();

  graphGrid->SetCells(VTK_LINE,meshLines);
  meshLines->Delete();
}

/**
 * @brief Program main.
 * @param argc argument counter
 * @param argv argument vector
 * @return rc the return code
 */
int main(int argc, char **argv)
{

  // STEP 0: Read table data
  if( argc != 2 )
    {
    std::cerr << "USAGE: mtreeviz <mtree.dat>\n";
    return EXIT_FAILURE;
    }

  std::string fileName( argv[1] );
  if( !ReadData( fileName ) )
    {
    std::cerr << "ERROR: could not read data successfully!\n";
    return EXIT_FAILURE;
    }

  // STEP 1: Build graph
  vtkSmartPointer<vtkMutableDirectedGraph> g =
   vtkSmartPointer<vtkMutableDirectedGraph>::New();
  BuildGraph( g );

  // STEP 3: Write poly data
  vtkUnstructuredGridWriter *writer = vtkUnstructuredGridWriter::New();
  writer->SetFileName( "ugrid.vtk" );
  writer->SetInputData( graphGrid );
  writer->Write();
  writer->Delete();

  // STEP 3: Visualize the graph
  vtkSmartPointer<vtkGraphLayoutView> graphLayoutView =
    vtkSmartPointer<vtkGraphLayoutView>::New();

  graphLayoutView->SetVertexScalarBarVisibility(true);
  graphLayoutView->AddRepresentationFromInput(g);
  graphLayoutView->SetColorVertices(true);
  graphLayoutView->SetVertexColorArrayName("HaloDensity");
  graphLayoutView->ResetCamera();
  graphLayoutView->Render();
  graphLayoutView->GetInteractor()->Start();

  graphGrid->Delete();
  return EXIT_SUCCESS;
}
