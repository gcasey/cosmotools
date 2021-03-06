cmake_minimum_required(VERSION 2.8)
project(pv-plugins)

## Find ParaView build and import cmake configurations, macros, etc.
find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})
include(ParaViewBranding)
include(ParaViewPlugins)

## Set the plugin name and version
set(PLUGIN_NAME PVCosmologyToolsPlugins)
set(PLUGIN_VERSION "1.0")

## Setup include directory paths
include_directories(
    ${PROJECT_SOURCE_DIR}/readers
    ${PROJECT_SOURCE_DIR}/filters
    ${CMAKE_SOURCE_DIR}/common
    ${CMAKE_SOURCE_DIR}/algorithms/sfprobe
    ${CMAKE_SOURCE_DIR}/algorithms/tess/lib
    ${PARAVIEW_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
    )

## Note we need to make sure that we put this include directory last so
## that the CosmologyTools halo-finder will be used and not the old out-dated
## version of the code that is available in ParaView.
include_directories(BEFORE ${PACKAGE_halofinder_INCLUDE_DIRS})
include_directories(BEFORE ${PACKAGE_tess_INCLUDE_DIRS})
include_directories(BEFORE ${PACKAGE_sfprobe_INCLUDE_DIRS})

## Set the server-manager XML
set(SMXML
    resources/Cosmo64Reader.xml
    resources/DensityProfile.xml
    resources/LANLHaloFinder.xml
    resources/StructureFormationProbe.xml
    resources/VoronoiFilter.xml
    resources/vtkPVoronoiReader.xml
    )

## Set the server-manager sources
set(SMSRC
    filters/VoronoiFilter.cxx
    filters/vtkCosmoDensityProfile.cxx
    filters/vtkMinkowskiFilter.cxx
    filters/vtkPLANLHaloFinder.cxx
    filters/vtkPMergeConnected.cxx
    filters/vtkPStructureFormationProbe.cxx
    readers/vtkPCosmo64Reader.cxx
    )

if(${ENABLE_QHULL})
set(SMSRC
    ${SMSRC}
    readers/vtkPVoronoiReader.cxx)
endif()

## Set the gui-resource_files
set(GUIXML
    resources/CosmologyToolsFilters.xml
    resources/CosmologyToolsReaders.xml
    )

## Add halo-finders to the project if that plug-in is built
set(RequiredLibs
    ${PACKAGE_halofinder_LIBRARIES}
    ${PACKAGE_tess_LIBRARIES}
    ${PACKAGE_sfprobe_LIBRARIES}
    ${MPI_C_LIBRARIES}
    ${VTK_LIBRARIES}
    )

## Add the paraview plugin
add_paraview_plugin(${PLUGIN_NAME} ${PLUGIN_VERSION}
    GUI_RESOURCE_FILES ${GUIXML}
    SERVER_MANAGER_XML ${SMXML}
    SERVER_MANAGER_SOURCES ${SMSRC}
    )
target_link_libraries(${PLUGIN_NAME} ${RequiredLibs})
