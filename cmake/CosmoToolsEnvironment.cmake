## Set the build type
set(CMAKE_BUILD_TYPE Debug CACHE STRING
  "one of: Release, Debug, RelWithDebInfo or MinSizeRel")

## Set compiler flags, etc.
include(CompilerExtras)

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

## This is needed for UINT64_C in GenericIO
add_definitions(-D__STDC_CONSTANT_MACROS)

## Choose static or shared libraries.
option(BUILD_SHARED_LIBS "Build shared libraries." OFF)
if(NOT BUILD_SHARED_LIBS)
 ## If we are building statically, make everything PICable so that we can
 ## link with shared libraries
 if(CMAKE_COMPILER_IS_GNUCC)
   set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} -fPIC)
   set(CMAKE_C_FLAGS ${CMAKE_C_FLAGS} -fPIC)
 endif() # END if compiler is GNU
endif() # END if not building shared libraries

## Enable framework statistics
option(ENABLE_FRAMEWORK_STATISTICS "Enable Framework Statistics" OFF)
if(${ENABLE_FRAMEWORK_STATISTICS})
  add_definitions(-DENABLESTATS)
endif()

option(ENABLE_THIRDPARTY_SQLITE "Enable SQlite" OFF)
if(${ENABLE_THIRDPARTY_SQLITE})
 add_definitions(-DSQLITE)
endif()
