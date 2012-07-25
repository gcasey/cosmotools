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
#include "vtkDataSetAttributes.h"
#include "vtkDoubleArray.h"
#include "vtkGraphLayout.h"
#include "vtkGraphLayoutView.h"
#include "vtkGraphToPolyData.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkTreeLayoutStrategy.h"

struct TreeNodeInfo {

  // Data read from file
  float HaloNumber;
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
    pos[0] = static_cast<double>(HaloNumber);
    pos[1] = this->TimeStep;
    pos[2] = 0.0;
    }

  /**
   * @return HashCode for
   */
 std::string GetHashCode()
   {
   std::ostringstream oss;
   oss << this->HaloNumber << "." << this->TimeStep;
   return( oss.str() );
   }

 /**
  * @brief Prints the info to the console.
  * @note For debugging.
  */
 void Println()
   {
   std::cout << "=======\n";
   std::cout << "HALO ID: "    << this->HaloNumber       << std::endl;
   std::cout << "DESCENDANT: " << this->DescendantNumber << std::endl;
   std::cout << "Time-step: "  << this->TimeStep         << std::endl;
   std::cout << "Volume: "     << this->Volume           << std::endl;
   std::cout << "PosX: "       << this->PosX             << std::endl;
   std::cout << "PosY: "       << this->PosY             << std::endl;
   std::cout << "PosZ: "       << this->PosZ             << std::endl;
   std::cout << "VelX: "       << this->VelX             << std::endl;
   std::cout << "VelY: "       << this->VelY             << std::endl;
   std::cout << "VelZ: "       << this->VelZ             << std::endl;
   std::cout << "hashcode: "   << this->GetHashCode()    << std::endl;
   std::cout.flush();
   }

};

// Stores all the nodes in the graph
std::vector< TreeNodeInfo > GraphNodes;

// Mapping of tree nodes to graph
std::map< std::string, int > TreeNode2Graph;

/**
 * @brief
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

    std::string hashCode = tnode.GetHashCode();
    if( TreeNode2Graph.find(hashCode) != TreeNode2Graph.end() )
      {
      std::cerr << "ERROR: Duplicate Entries in TreeNode2Graph\n";
      return false;
      }

    TreeNode2Graph[ hashCode ] = GraphNodes.size()-1;
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

  vtkDoubleArray *volumes = vtkDoubleArray::New();
  volumes->SetName("HaloDensity");
  volumes->SetNumberOfComponents( 1 );
  volumes->SetNumberOfTuples( GraphNodes.size() );

  double *vptr = static_cast<double*>(volumes->GetPointer(0));

  // STEP 0: Insert nodes
  double pos[3];
  for( int i=0; i < GraphNodes.size(); ++i )
    {
    GraphNodes[i].GetGraphNodePosition(pos);
    vptr[ i ] = GraphNodes[ i ].Volume;

    vtkIdType idx = g->AddVertex();
//    g->GetPoints()->SetPoint( idx, pos);
    } // END for all graph nodes
  g->GetVertexData()->AddArray( volumes );
  volumes->Delete();

  // STEP 1: Create edges
  for( int i=0; i < GraphNodes.size(); ++i )
    {
    if( static_cast<int>(GraphNodes[i].DescendantNumber) != -1)
      {
      g->AddEdge( i, GraphNodes[i].DescendantNumber);
      }
    } // END for all graph nodes
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

  // STEP 2: Visualize the graph
  vtkSmartPointer<vtkGraphLayoutView> graphLayoutView =
    vtkSmartPointer<vtkGraphLayoutView>::New();

  // TODO: Ask Jeff about how to use the Tree Layout strategy -- The problem
  // is that Boost is required but, VTK_USE_BOOST is not available as an option
  // that can be turned on. When manually getting rid of the #ifdef VTK does
  // not compile.
  //
//  graphLayoutView->SetLayoutStrategyToCosmicTree();
  //


  graphLayoutView->AddRepresentationFromInput(g);
  graphLayoutView->SetColorVertices(true);
  graphLayoutView->SetVertexColorArrayName("HaloDensity");
  graphLayoutView->ResetCamera();
  graphLayoutView->Render();
  graphLayoutView->GetInteractor()->Start();

  return EXIT_SUCCESS;
}
