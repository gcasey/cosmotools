## Set the minimum version for cmake to 2.8
cmake_minimum_required(VERSION 2.8.11)

## Set the project name and enable C/C++ support by default
project(CosmoTools C CXX)

## Set path for cmake modules
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

include(CosmoToolsEnvironment)
include(CosmoToolsMacros)
include(CosmoToolsDependencies)
include(CosmoToolsTypes)

include_directories(${HEADER_INCLUDES_DIRECTORY})

## list of required libraries for the entire package
set(CosmoToolsRequiredLibs ${MPI_C_LIBRARIES} ${GENERIC_IO_LIBRARIES})

## Set the list of module that will be built. Note the order here matters!
set(modules
    common
    algorithms
    framework
    )

## Builds C [and Fortran] interface
if(BUILD_SIMULATION_INTERFACE)
  list(APPEND modules interfaces)
endif()

## Loop through all packages
foreach(mod ${modules})
  add_subdirectory(${mod})
endforeach() ## END for each module

## Build a single cosmotools library
if(BUILD_SINGLE_LIBRARY)
    if(BUILD_SHARED_LIBS)
        add_library(cosmotools SHARED
            $<TARGET_OBJECTS:common_obj>
            $<TARGET_OBJECTS:halofinder_obj>
            $<TARGET_OBJECTS:framework_obj>
            )
    else()
        add_library(cosmotools STATIC
            $<TARGET_OBJECTS:common_obj>
            $<TARGET_OBJECTS:halofinder_obj>
            $<TARGET_OBJECTS:framework_obj>
            )
    endif() # END if Build shared libraries
    target_link_libraries(cosmotools ${CosmoToolsRequiredLibs})
endif() # END if Build single library

option(BUILD_PV_PLUGINS "Build paraview plugins" OFF)

if(BUILD_PV_PLUGINS)
  add_subdirectory(pv-plugins)
endif()
