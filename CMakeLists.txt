## Set the minimum version for cmake to 2.8
cmake_minimum_required(VERSION 2.8)

## Set the project name and enable C/C++ and Fortran support
project(CosmologyTools C CXX Fortran)

## Set path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(CosmologyToolsEnvironment)
include(CosmologyToolsTypes)
include(GetPackageSources)

## Set the list of all the cosmotools sources
set(cosmosrclist)

## Set the list of cosmology packages to be build. Note the order here matters!
set(cosmopackages
    common
    sfprobe
    halofinder
    tess
    mergertree
    )

## Sets whether or not the simulation interface will be built
option(BUILD_SIMULATION_INTERFACE "Build Simulation inteface" OFF)
if(${BUILD_SIMULATION_INTERFACE})
  set(cosmopackages
    ${cosmopackages}
    framework
    interfaces
    )
endif()

## Sets whether or not ParaView plugins will be built
option(BUILD_PV_PLUGINS "Build ParaView Plugins" OFF)
if(${BUILD_PV_PLUGINS})
 set(cosmopackages
     ${cosmopackages}
     pv-plugins
     )
endif()

## Sets whether or not VTK programs will be built
option(BUILD_VTK_PROGRAMS "Build VTK Programs" OFF)
if(${BUILD_VTK_PROGRAMS})
 find_package(VTK COMPONENTS
    vtkCommonDataModel
    vtkFiltersCosmo
    vtkIOXML
    )
 set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${VTK_DEFINITIONS})
 include_directories(SYSTEM ${VTK_INCLUDE_DIRS})
endif()


## Loop through all packages
set(emptyList)
foreach(pkg ${cosmopackages})
  ## These are declared here and populated by each package
  ## NOTE: These variables must be set before the call to add_subdirectory
  set(PACKAGE_${pkg}_INCLUDE_DIRS ${emptyList}
      CACHE INTERNAL "package include directories" FORCE)
  set(PACKAGE_${pkg}_LIBRARIES ${emptyList}
      CACHE INTERNAL "package libraries" FORCE)
  set(PACKAGE_${pgk}_SOURCES ${emptyList}
      CACHE INTERNAL "package sources" FORCE)

  ## add package to build
  add_subdirectory(${pkg})

  ## append package sources to the cosmosrclist
  list(APPEND cosmosrclist ${PACKAGE_${pkg}_SOURCES})
endforeach()

## Build single library
if(BUILD_SHARED_LIBS)
 add_library(cosmotools SHARED ${cosmosrclist})
else()
 add_library(cosmotools STATIC ${cosmosrclist})
endif()
target_link_libraries(cosmotools ${MPI_LIBRARIES})
