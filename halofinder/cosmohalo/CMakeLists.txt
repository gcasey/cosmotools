cmake_minimum_required(VERSION 2.8)
project(cosmohalo)

## Specicy global include dirs and global libraries for this project
## such that, they can be used within other projects.
set(PACKAGE_halofinder_INCLUDE_DIRS
    ${PACKAGE_halofinder_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "package include directories")
set(PACKAGE_halofinder_LIBRARIES
   ${PACKAGE_halofinder_LIBRARIES}
   cosmohalo
   CACHE INTERNAL "package libraries")

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Choose whether to use int64_t or int32_t
option(ID64 "Use int64_t for particles and halo ids" ON)
if( ${ID64} )
 add_definitions(-DID_64)
endif()

## Choose whether to use double or float for position,velocity and potential
option(POSVEL64 "Use double precission for position,velocity and potential" ON)
if( ${POSVEL64} )
 add_definitions(-DPOSVEL_64)
endif()

## Choose whether to use double or float for grid types
option(GRID64 "Use double precision for grid types" ON)
if( ${GRID64} )
 add_definitions(-DGRID_64)
endif()

set(Cosmo_LIBS  ${MPI_LIBRARIES})

## List of executables
set(Cosmo_PROGRAMS
  BinHalos.cxx
  ForceRandom.cxx
  ForceSpeed.cxx
  NToN.cxx
  OneToN.cxx
  )

## Set cosmo sources
set(Cosmo_SOURCES
  BHForceTree.cxx
  BHTree.cxx
  bigchunk.cxx
  ChainingMesh.cxx
  CosmoHaloFinder.cxx
  CosmoHaloFinderP.cxx
  dims.cxx
  FOFHaloProperties.cxx
  ForceLaw.cxx
  ForceTree.cxx
  GridExchange.cxx
  HaloCenterFinder.cxx
  HaloFinderInput.cxx
  InitialExchange.cxx
  Message.cxx
  Partition.cxx
  ParticleDistribute.cxx
  ParticleExchange.cxx
  RCBForceTree.cxx
  RCOForceTree.cxx
  SODHalo.cxx
  SubHaloFinder.cxx
  Timer.cxx
  Timings.cxx
  )

## Create the libraries
if(BUILD_SHARED_LIBS)
 add_library(cosmohalo SHARED ${Cosmo_SOURCES})
else()
 add_library(cosmohalo STATIC ${Cosmo_SOURCES})
endif()
target_link_libraries(cosmohalo ${MPI_LIBRARIES})

## Compile each program
foreach(prgm ${Cosmo_PROGRAMS})
  get_filename_component(appName ${prgm} NAME_WE)
  add_executable( ${appName} ${prgm} )
  target_link_libraries(${appName} ${Cosmo_LIBS} cosmohalo ${MPI_LIBRARIES})
endforeach()
