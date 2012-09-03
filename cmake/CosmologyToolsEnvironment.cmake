## Set the build type
set(CMAKE_BUILD_TYPE Debug CACHE STRING
  "one of: Release, Debug, RelWithDebInfo or MinSizeRel" FORCE)

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
if(NOT BUILD_SHARED_LIBS)
 ## If we are building statically, make everything PICable so that we can
 ## link with shared libraries
 add_definitions(-fPIC)
endif()

## Choose whether to build tests or not
option(BUILD_TESTING "Build Tests." OFF)
