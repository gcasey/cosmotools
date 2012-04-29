## Set the minimum version for cmake to 2.8
cmake_minimum_required(VERSION 2.8)

## Set the project name and enable C/C++ and Fortran support
project(HACCTools C CXX Fortran)

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


## Loop thrhough all package directories and provide an option for which
## package to be enabled, i.e., set to ON. All enabled packages will be built.
set(subdirs)
file(GLOB pkgs *)
foreach(pkgpath ${pkgs})
 if(IS_DIRECTORY ${pkgpath} )
  ## Extract relative path to directory
  string(REPLACE "/" ";" pkglist ${pkgpath})
  list(REVERSE pkglist)
  list(GET pkglist 0 pkg)

  string(COMPARE NOTEQUAL ${pkg} "common" isOptionalPackage)
  if( ${isOptionalPackage} )
    option(BUILD_${pkg} "Build ${pkg}." OFF)
    if(BUILD_${pkg})
        list(APPEND subdirs ${pkg})
    endif() # end BUILD_${pkg}
  endif() # end is optional package
 endif() # end is_directory
endforeach(pkgpath)

add_subdirectory(common)

## Add each package
foreach(subdir ${subdirs})
  add_subdirectory(${subdir})
endforeach()


