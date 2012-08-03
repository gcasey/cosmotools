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

if( ${USECPP0X} )
 if( ${USEINT64} )
  add_definitions(-DID_64)
 endif()
endif()

if(${DOUBLE_PRECISION})
  add_definitions(-DGRID_64 -DPOSVEL_64)
endif()

set(Cosmo_LIBS  ${MPI_LIBRARIES})

## List of executables
set(Cosmo_PROGRAMS
  BinHalos.cxx
  NToN.cxx
  OneToN.cxx
  )

## Set cosmo sources
set(Cosmo_SOURCES
  BHTree.cxx
  bigchunk.cxx
  ChainingMesh.cxx
  CosmoHaloFinder.cxx
  CosmoHaloFinderP.cxx
  dims.cxx
  FOFHaloProperties.cxx
  GridExchange.cxx
  HaloCenterFinder.cxx
  HaloFinderInput.cxx
  InitialExchange.cxx
  Message.cxx
  Partition.cxx
  ParticleDistribute.cxx
  ParticleExchange.cxx
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
