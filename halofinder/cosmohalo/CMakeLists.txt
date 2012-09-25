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

set(Cosmo_LIBS  ${MPI_C_LIBRARIES})

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

## Specify global source list for this project
#set(PACKAGE_halofinder_SOURCES
#    ${PACKAGE_halofinder_SOURCES}
#    ${Cosmo_SOURCES}
#    CACHE INTERNAL "package sources")
GetPackageSources(halofinder ${PROJECT_SOURCE_DIR} "${Cosmo_SOURCES}")

## Create the libraries
if(BUILD_SHARED_LIBS)
 add_library(cosmohalo SHARED ${Cosmo_SOURCES})
else()
 add_library(cosmohalo STATIC ${Cosmo_SOURCES})
endif()
target_link_libraries(cosmohalo ${MPI_C_LIBRARIES})

## Compile each program
foreach(prgm ${Cosmo_PROGRAMS})
  get_filename_component(appName ${prgm} NAME_WE)
  add_executable( ${appName} ${prgm} )
  target_link_libraries(${appName} ${Cosmo_LIBS} cosmohalo ${MPI_C_LIBRARIES})
endforeach()
