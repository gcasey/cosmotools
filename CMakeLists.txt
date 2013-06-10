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
set(CosmoToolsRequiredLibs ${MPI_C_LIBRARIES})
if(${ENABLE_DIY})
 list(APPEND CosmoToolsRequiredLibs ${DIY_LIBRRARIES})
endif()


## Set the list of cosmology packages to be build. Note the order here matters!
set(modules
    common
    )

## Sets whether or not the simulation interface will be built
#option(BUILD_SIMULATION_INTERFACE "Build Simulation inteface" OFF)
#if(${BUILD_SIMULATION_INTERFACE})
#  set(cosmopackages
#    ${cosmopackages}
#    framework
#    interfaces
#    )
#endif()

## Sets whether or not ParaView plugins will be built
#option(BUILD_PV_PLUGINS "Build ParaView Plugins" OFF)
#if(${BUILD_PV_PLUGINS})
# set(cosmopackages
#     ${cosmopackages}
#     pv-plugins
#     )
#endif()

## Loop through all packages
set(emptyList)
foreach(module_m ${modules})
  add_subdirectory(${module_m})
endforeach()

if(${BUILD_SINGLE_LIBRARY})
    if(BUILD_SHARED_LIBS)
        add_library(cosmotools SHARED $<TARGET_OBJECTS:common_obj>)
    else()
        add_library(cosmotools STATIC $<TARGET_OBJECTS:common_obj>)
    endif()
    target_link_libraries(cosmotools ${CosmoToolsRequiredLibs})
endif()

#if(${BUILD_SINGLE_LIBRARY})
# if(BUILD_SHARED_LIBS)
#  add_library(cosmotools SHARED ${cosmosrclist})
# else()
#  add_library(cosmotools STATIC ${cosmosrclist})
# endif()
# target_link_libraries(cosmotools ${CosmoRequiredLibs})
#endif()
