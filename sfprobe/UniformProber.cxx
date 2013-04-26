#include "UniformProber.h"

#include "ExtentUtilities.h"
#include "StructureFormationProbe.h"
#include "VirtualGrid.h"

#ifdef USEDAX
 // Dax includes
#  include "dax/ConstArray.h"
#  include "dax/MapPointsToGrid.h"
#  include "dax/MapTetsToPoints.h"
#  include "dax/TetStructure.h"

#  include <dax/CellTag.h>
#  include <dax/cont/arg/ExecutionObject.h>
#  include <dax/cont/DeviceAdapter.h>
#  include <dax/cont/Scheduler.h>
#  include <dax/cont/Timer.h>
#  include <dax/Extent.h>
#  include <dax/math/Compare.h>

#  include <vector>
#endif

#include <iostream>


#ifdef USEDAX
namespace detail
{
template< class DeviceAdapterTag, class ArrayValueType, class ArrayContainerType>
void InvokeMapPointsToGrid(
      const dax::cont::Scheduler<DeviceAdapterTag>& scheduler,
      dax::cont::UniformGrid<DeviceAdapterTag> probeGrid,
      cosmologytools::VirtualGrid& vgrid,
      dax::cont::ArrayHandle< ArrayValueType, ArrayContainerType, DeviceAdapterTag> vgridCellId)
  {
  //extract out the spacing, origin, and extents of virtual grid
  dax::Vector3 vgrid_spacing(vgrid.GetSpacing()[0],
                             vgrid.GetSpacing()[1],
                             vgrid.GetSpacing()[2]);
  const REAL* vgrid_bounds = vgrid.GetBounds();
  dax::Vector3 vgrid_origin( vgrid_bounds[0], vgrid_bounds[2], vgrid_bounds[4] );
  dax::Extent3 vgrid_ext( dax::make_Id3(0,0,0),
                          dax::Id3(vgrid.GetDimensions()[0] - 1,
                                   vgrid.GetDimensions()[1] - 1,
                                   vgrid.GetDimensions()[2] - 1) );
  scheduler.Invoke( dax::worklet::MapPointsToGrid(),
                    probeGrid.GetPointCoordinates(),
                    vgrid_origin,
                    vgrid_spacing,
                    vgrid_ext,
                    vgridCellId);

  //manually sort the points
  dax::math::SortLess comparisonFunctor;
  dax::cont::internal::DeviceAdapterAlgorithm<DeviceAdapterTag>::Sort(
                                              vgridCellId, comparisonFunctor);

  }
}
#endif

namespace cosmologytools {

UniformProber::UniformProber( REAL origin[3], REAL spacing[3], INTEGER ext[6] )
{
  this->Origin[0] = origin[0];
  this->Origin[1] = origin[1];
  this->Origin[2] = origin[2];

  this->Spacing[0] = spacing[0];
  this->Spacing[1] = spacing[1];
  this->Spacing[2] = spacing[2];

  this->Extents[0] = ext[0];
  this->Extents[1] = ext[1];
  this->Extents[2] = ext[2];
  this->Extents[3] = ext[3];
  this->Extents[4] = ext[4];
  this->Extents[5] = ext[5];

  //compute the number of points in the grid
  this->NumberOfPoints =
    cosmologytools::ExtentUtilities::ComputeNumberOfNodes(ext);

  this->NumberOfStreams = new INTEGER[this->NumberOfPoints];
  this->Rho = new REAL[this->NumberOfPoints];
}

UniformProber::~UniformProber()
{
  if(this->NumberOfStreams)
    {
    delete[] this->NumberOfStreams;
    this->NumberOfStreams  = NULL;
    }

  if(this->Rho)
    {
    delete[] this->Rho;
    this->Rho = NULL;
    }
}

void UniformProber::ComputePoint(INTEGER index, REAL xyz[3] )
{
  INTEGER ijk[3];
  cosmologytools::ExtentUtilities::GetStructuredCoordinates(index,
                                                            this->Extents,
                                                            ijk);
  xyz[0] = this->Origin[0] + this->Spacing[0]* ijk[0];
  xyz[0] = this->Origin[1] + this->Spacing[1]* ijk[1];
  xyz[0] = this->Origin[2] + this->Spacing[2]* ijk[2];
}

void UniformProber::RunProber( cosmologytools::StructureFormationProbe * probe )
{
#ifndef USEDAX
  this->RunSerialProber(probe);
#else
  this->RunDaxProber(probe);
#endif
}

void UniformProber::RunSerialProber( cosmologytools::StructureFormationProbe * probe )
{
  #pragma omp parallel for
  for(int pntIdx=0; pntIdx < this->NumberOfPoints; ++pntIdx )
    {
    REAL pnt[3];
    this->ComputePoint(pntIdx,pnt);

    INTEGER nStream = 0;
    REAL lrho = 0.0;
    probe->ProbePoint(pnt,nStream,lrho);

    this->NumberOfStreams[ pntIdx ] = nStream;
    this->Rho[ pntIdx ] = lrho;
    } // END for all points

  //Statistics:
  std::cout << "NumPoints probed=" << probe->GetNumPointsProbed();
  std::cout << std::endl;
  std::cout << "NumTets checked=" << probe->GetNumTetsChecked();
  std::cout << std::endl;
  std::cout << "Avg. NumTets per point=";
  std::cout << probe->GetNumTetsChecked()/probe->GetNumPointsProbed();
  std::cout << std::endl;
}

void UniformProber::RunDaxProber( cosmologytools::StructureFormationProbe * probe )
{
#ifdef USEDAX
  dax::cont::Scheduler<> scheduler;
  dax::cont::Timer<> timer;

  //STEP 1:
  // construct the virtual grid, which will classify each tet to a bucket
  //We set the dims to be 8 in each direction as a guestimate on a reasonable
  //number to decompose them problem at. We want it to be fairly coarse
  //for TBB, unkown what performs best for Cuda.
  INTEGER vgrid_dims[3] = {8,8,8};
  cosmologytools::VirtualGrid vgrid;
  vgrid.SetDimensions(vgrid_dims);
  vgrid.RegisterMesh( probe->GetEulerMesh() );


  //STEP 2:
  // Build a Dax uniform grid the same size as our self that we can
  // use to as the probe points

  //map each point in the probe set to the virtual grid
  dax::cont::UniformGrid<> probeGrid;
  probeGrid.SetOrigin( dax::make_Vector3(
                          this->Origin[0],this->Origin[1],this->Origin[2]) );
  probeGrid.SetSpacing( dax::make_Vector3(
                          this->Spacing[0],this->Spacing[1],this->Spacing[2]) );
  probeGrid.SetExtent(  dax::make_Id3(this->Extents[0], this->Extents[2], this->Extents[4]),
                        dax::make_Id3(this->Extents[1], this->Extents[3], this->Extents[5]));


  //STEP 3:
  //Ask Dax to map each of the points in the created uniform grid, to a bucket
  //in the VirtualGrid. The vgridCellId holds the result, which is a Tuple
  //where the first element is the virtual grid bucket and the second element
  //is the point id. This also sorts the Tuples so that all the points
  //that have the same bucket are continuous in memory
  typedef dax::cont::ArrayHandle< dax::Tuple<dax::Id, 2> >  TupleHandleType;

  TupleHandleType vgridCellId;
  detail::InvokeMapPointsToGrid(scheduler, probeGrid, vgrid, vgridCellId);

  //STEP 4: Get the sorted results back from Dax.
  //pull down the mapping
  typedef TupleHandleType::PortalConstControl PortalType;
  PortalType values = vgridCellId.GetPortalConstControl();


  //STEP 5: Zero out the rho and num streams arrays
  std::fill(this->NumberOfStreams,this->NumberOfStreams+NumberOfPoints,0);
  std::fill(this->Rho,this->Rho+NumberOfPoints,0);


  //STEP 6: We construct the Tet Search structure for the final worklet
  // WE Extract all the tets from the Euler mesh to do this
  //these will be the tets we will search
  std::vector<dax::Scalar> nodes;
  std::vector<dax::Id> tetConn;
  std::vector<dax::Scalar> volumes;
  probe->GetEulerMesh(nodes,tetConn,volumes);

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

  dax::exec::TetStructure tetStructure(tetGrid,volumeHandle);

  //STEP 7:
  //Remove all cells that have been mapped to a negative cell id
  //which means they lay outside the virtual grid
  const int maxBucketNum = vgrid.GetNumberOfBuckets();
  dax::Id min_valid_id=0;
  while( values.Get(min_valid_id++)[0] < 0 );

  //STEP 8:
  //Iterate over each bucket getting collecting all the points that
  //are inside of it. Dispatch a dax worklet that calculates out for each
  //point the rho and num streams

  for(dax::Id i=min_valid_id; i < this->NumberOfPoints-1; ++i)
    {
    //STEP 8 A:
    //make vector to hold the point coordinates and point indices
    //for each point in this bucket
    std::vector< dax::Vector3 > subCoords;
    std::vector< dax::Id > writeIndices;
    subCoords.reserve(128); //random guess
    writeIndices.reserve(128); //random guess

    //STEP 8 B:
    //Continue looping over points add them to the vectors while they share
    //the same cellIndex
    const dax::Id cellIndex = values.Get(i)[0];
    const dax::Id pointIndex = values.Get(i)[1];
    subCoords.push_back( probeGrid.ComputePointCoordinates( pointIndex  ) );
    writeIndices.push_back( pointIndex );
    while( i < this->NumberOfPoints-1 && cellIndex == values.Get(i+1)[0])
      {
      const dax::Id nextPointIndex = values.Get(i+1)[1];
      subCoords.push_back( probeGrid.ComputePointCoordinates(nextPointIndex) );
      writeIndices.push_back( nextPointIndex );
      ++i;
      }

    //if we only have a single point just use the serial code
    std::size_t numPointsInBucket = subCoords.size();
    if(numPointsInBucket == 1)
      {
      REAL xyz[3] = {subCoords[0][0],subCoords[0][1],subCoords[0][2]};
      INTEGER nStream = 0;
      REAL lrho = 0.0;
      probe->ProbePoint(xyz,nStream,lrho);

      this->NumberOfStreams[ writeIndices[0] ] = static_cast<int>(nStream);
      this->Rho[ writeIndices[0] ]     = static_cast<double>(lrho);
      } // END for all points
    else if (numPointsInBucket > 1)
      {
      //STEP 8 C:
      //Finally we have the information needed to dispatch to Dax
      //so invoke on the subset of points and tets that exist in the bucket

      //get all the tets in the bucket, bucket id will always be valid
      std::set<dax::Id> const& tetSet =  vgrid.GetBucket( cellIndex );
      std::vector<dax::Id> tetIds(tetSet.begin(),tetSet.end());

      //wrap coords in an array handle
      dax::cont::ArrayHandle< dax::Vector3 > subCoordsHandle =
          dax::cont::make_ArrayHandle( subCoords );

      //wrap the tetIds in the bucket in a exec object
      dax::cont::ArrayHandle< dax::Id > tetIdHandle =
          dax::cont::make_ArrayHandle( tetIds );
      dax::exec::ConstArray<dax::Id> tetIdsToSearch(tetIdHandle);

      //setup handles to the output data
      dax::cont::ArrayHandle<dax::Id> streamsHandle;
      dax::cont::ArrayHandle<dax::Scalar> rhoHandle;

      scheduler.Invoke(dax::worklet::MapTetsToPoints(),
          subCoordsHandle,
          tetStructure,
          tetIdsToSearch,
          streamsHandle,
          rhoHandle);

      //STEP 8 D:
      //Write the resulting rho and num streams back to the host vectors
      typedef std::vector<dax::Id>::const_iterator iterator;
      const dax::Id* streamIt = streamsHandle.GetPortalConstControl().GetIteratorBegin();
      const dax::Scalar* rhoIt = rhoHandle.GetPortalConstControl().GetIteratorBegin();
      for(iterator j = writeIndices.begin(); j != writeIndices.end(); ++j )
        {
        this->NumberOfStreams[*j] = *streamIt;
        this->Rho[*j] = static_cast<double>(*rhoIt);
        ++streamIt;
        ++rhoIt;
        }
      }
    }
#endif
}

}
