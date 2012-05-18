## Set the minimum version for cmake to 2.8
cmake_minimum_required(VERSION 2.8)

## Set the project name and enable C/C++ and Fortran support
project(CosmologyTools C CXX Fortran)

## Set path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

## Set the build type
set(CMAKE_BUILD_TYPE "Debug" CACHE STRING
  "one of: Release, Debug, RelWithDebInfo or MinSizeRel")

## Set the path where all the header will be stored
set(HEADER_INCLUDES_DIRECTORY
    ${PROJECT_BINARY_DIR}/include
    CACHE PATH
    "Directory where all headers will live"
    )

## Set the path where all the libraries will be stored
set(LIBRARY_OUTPUT_PATH
    ${PROJECT_BINARY_DIR}/libs
    CACHE PATH
    "Directory where compiled libraries will live"
    )

## Set the path where all the executables will go
set(EXECUTABLE_OUTPUT_PATH
    ${PROJECT_BINARY_DIR}/bin
    CACHE PATH
    "Directory where executables will live"
    )

## Set the Fortran module directory
set(CMAKE_Fortran_MODULE_DIRECTORY
    ${PROJECT_BINARY_DIR}/FortranModules
    CACHE PATH
    "Directory where all Fortran modules will live"
    )

## Mark as advanced
mark_as_advanced(
    LIBRARY_OUTPUT_PATH
    EXECUTABLE_OUTPUT_PATH
    CMAKE_Fortran_MODULE_DIRECTORY)

## Choose whether to turn on verbosity or not
option(MAKEFILE_VERBOSE "Turn on verbose makefiles." OFF)
if(${MAKEFILE_VERBOSE})
 set(CMAKE_VERBOSE_MAKEFILE ON)
endif()

## MPI is a required dependency
include(FindMPI REQUIRED)
include_directories(
    common
    ${MPI_INCLUDE_PATH})
add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

## Choose static or shared libraries.
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)

## Choose whether to use double-precision or not
option(DOUBLE_PRECISION "Enable double-precision." OFF)
if(${DOUBLE_PRECISION})
 add_definitions(-DUSE_DOUBLE_PRECISION)
endif()

## Choose whether to use int64 or int32
option(USECPP0X "Use features from the upcoming ISO C++ standard" OFF)
if( ${USECPP0X} )
 add_definitions(-DUSECPP0X --std=c++0x)
 option(USEINT64 "Use int64 or int32" OFF)
 if( ${USEINT64} )
  add_definitions(-DUSEINT64 )
 endif()
endif()

## Set the list of cosmology packages to be build. Note the order here matters!
set(cosmopackages
    common
    sfprobe
    halofinder
    tess
    interfaces
    pv-plugins
    )

## Loop through all packages
set(emptyList)
foreach(pkg ${cosmopackages})
  ## These are declared here and populated by each package
  ## NOTE: These variables must be set before the call to add_subdirectory
  set(PACKAGE_${pkg}_INCLUDE_DIRS ${emptyList}
      CACHE INTERNAL "package include directories" FORCE)
  set(PACKAGE_${pkg}_LIBRARIES ${emptyList}
      CACHE INTERNAL "package libraries" FORCE)

  ## add package to build
  add_subdirectory(${pkg})

endforeach()
