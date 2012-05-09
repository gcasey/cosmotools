/*=========================================================================

 Program:   Visualization Toolkit
 Module:    vtkCosmoDensityProfile.h

 Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
 All rights reserved.
 See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

 This software is distributed WITHOUT ANY WARRANTY; without even
 the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
 PURPOSE.  See the above copyright notice for more information.

 =========================================================================*/
// .NAME vtkCosmoDensityProfile.h -- Interactive density profile tool
//
// .SECTION Description
// A concrete instance of vtkPolyDataAlgorithm that provides functionality for
// acquiring the density profile at a user-prescribed location. The user
// provides a center, a radius and N, the number of points, p_i, to subdivide
// the vector from the user-supplied center with a magnitude equal to the
// radius. Con-centric spheres are created through each point p_i. For each
// sphere the number of particles within the sphere is computed and stored in
// a vtkTable instance. The output of this filter consists of the set of spheres
// and the table. The input to the filter is a vtkUnstructuredGrid which is
// used to represent particle datasets.
//
// .SECTION See Also
//  vtkPolyDataAlgorithm

#ifndef VTKCOSMODENSITYPROFILE_H_
#define VTKCOSMODENSITYPROFILE_H_

#include "vtkMultiBlockDataSetAlgorithm.h"
#include <vector> // For STL vector
#include <list> // For STL list

// Forward declarations
class vtkAlgorithm;
class vtkAlgorithmOutput;
class vtkIndent;
class vtkMultiProcessController;
class vtkInformation;
class vtkInformationVector;
class vtkUnstructuredGrid;
class vtkSphere;

class VTK_EXPORT vtkCosmoDensityProfile : public vtkMultiBlockDataSetAlgorithm
{
public:
  static vtkCosmoDensityProfile* New();
  vtkTypeMacro(vtkCosmoDensityProfile,vtkMultiBlockDataSetAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get radius
  vtkSetMacro(Radius,int);
  vtkGetMacro(Radius,int);

  // Description:
  // Set/Get the number of points.
  vtkSetMacro(NumberOfSpheres,int);
  vtkGetMacro(NumberOfSpheres,int);

  // Description:
  // Set/Get the center
  vtkSetVector3Macro(Center,double);
  vtkGetVector3Macro(Center,double);

  // Description:
  // Set & Get the multi-process controller.
  vtkGetMacro(Controller,vtkMultiProcessController*);
  vtkSetMacro(Controller,vtkMultiProcessController*);

  // Description:
  // Method used by the pipeline to connect the FOF centers output from the
  // halo-finder to this algorithm.
  void SetFOFCentersConnection(vtkAlgorithmOutput* algOutput);

protected:
  vtkCosmoDensityProfile();
  virtual ~vtkCosmoDensityProfile();

  // Description:
  // Standard pipeline operations
  virtual int RequestData( vtkInformation *request,
      vtkInformationVector **inputVector,
      vtkInformationVector *outputVector );
  virtual int FillInputPortInformation(int port, vtkInformation *info);
  virtual int FillOutputPortInformation( int port, vtkInformation *info );

  // Description:
  // Generate con-centric spheres
  void GenerateSpheres(
       const int sphereIdx,
      double center[3],
      vtkMultiBlockDataSet *output);

  // Description:
  // Computes the density, i.e., the number of particles within each sphere
  // at the user-supplied region.
  void ComputeDensity(
      const int sphereIdx,
      double center[3],
      vtkUnstructuredGrid* particles);

  // Description:
  // Finds all the particles that are within the given sphere
  void FindParticlesInSphere(
      const int radiusIdx,
      double center[3],
      vtkUnstructuredGrid* particles,
      std::list< vtkIdType > &particleIds);

  // Description:
  // Computes the Halo FOF centers. Set the NumberOfSphereCenters equal to the
  // number of halos and populates the SphereCenters with the FOF centers array.
  void GetHaloFOFCenters(vtkUnstructuredGrid *fofCenters);

  // Description:
  // Computes a density profile at the given sphere center.
  void ProcessSphereCenter(
      int sphereIdx,
      vtkUnstructuredGrid *particles,
      vtkMultiBlockDataSet *spheres );

  int NumberOfSpheres; // The number of concentric spheres
  int Radius;          // The radius to use from the sphere center
  double Center[3];    // The sphere center to use when in interactive mode
  int NumberOfSphereCenters; // The number of sphere centers to process
  int UseFOFCenters; // A flag that indicates whether or not to process FOFs
  vtkMultiProcessController *Controller; // The controller used in parallel

//BTX
  std::vector<double> SphereCenters; // The list of sphere-centers to process
  std::vector<double> ConcentricRadii; // The radii of the concentric spheres
  std::vector<int> NumParticlesInSphere; // The density of each sphere.
//ETX

private:
  vtkCosmoDensityProfile(const vtkCosmoDensityProfile&); // Not implemented
  void operator=(const vtkCosmoDensityProfile&); // Not implemented
};

#endif /* VTKCOSMODENSITYPROFILE_H_ */
