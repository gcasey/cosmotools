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

#include "vtkPolyDataAlgorithm.h"

class vtkIndent;
class vtkMultiProcessController;
class vtkInformation;
class vtkInformationVector;

class VTK_EXPORT vtkCosmoDensityProfile : public vtkPolyDataAlgorithm
{
public:
  static vtkCosmoDensityProfile* New();
  vtkTypeMacro(vtkCosmoDensityProfile,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Set/Get radius
  vtkSetMacro(Radius,int);
  vtkGetMacro(Radius,int);

  // Description:
  // Set/Get the number of points.
  vtkSetMacro(NumberOfPoints,int);
  vtkGetMacro(NumberOfPoints,int);

  // Description:
  // Set/Get the center
  vtkSetVector3Macro(Center,double);
  vtkGetVector3Macro(Center,double);

  // Description:
  // Set & Get the multi-process controller.
  vtkGetMacro(Controller,vtkMultiProcessController*);
  vtkSetMacro(Controller,vtkMultiProcessController*);

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
  // Generate spheres
  void GenerateSpheres();

  int NumberOfPoints;
  int Radius;
  double Center[3];
  vtkMultiProcessController *Controller;

private:
  vtkCosmoDensityProfile(const vtkCosmoDensityProfile&); // Not implemented
  void operator=(const vtkCosmoDensityProfile&); // Not implemented
};

#endif /* VTKCOSMODENSITYPROFILE_H_ */
