## Set the build type
set(CMAKE_BUILD_TYPE Debug CACHE STRING
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
    ${HEADER_INCLUDES_DIRECTORY}
    common
    ${MPI_C_INCLUDE_PATH})
add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

## DIY is a required dependency, find it here
include(FindDIY REQUIRED)
include_directories(${DIY_INCLUDE_DIRS})


## Choose static or shared libraries.
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)
if(NOT BUILD_SHARED_LIBS)
 ## If we are building statically, make everything PICable so that we can
 ## link with shared libraries
 add_definitions(-fPIC)
endif()

option(ENABLE_FRAMEWORK_STATISTICS "Enable Framework Statistics" OFF)
if(${ENABLE_FRAMEWORK_STATISTICS})
  add_definitions(-DENABLESTATS)
endif()

option(ENABLE_DAX "Enable DAX Toolkit" OFF)
if(${ENABLE_DAX})
  add_definitions(-DUSEDAX)
  find_package(Dax REQUIRED)
  DaxConfigureTBB(REQUIRED)
  DaxConfigureCuda()
endif()

option(ENABLE_THIRDPARTY_SQLITE "Enable SQlite" OFF)
if(${ENABLE_THIRDPARTY_SQLITE})
 add_definitions(-DSQLITE)
endif()

## Choose whether to build tests or not
option(BUILD_TESTING "Build Tests." OFF)
