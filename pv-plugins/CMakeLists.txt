cmake_minimum_required(VERSION 2.8)
project(pv-plugins)

## Find ParaView build and import cmake configurations, macros, etc.
find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})

## Set the plugin name and version
set(PLUGIN_NAME PVCosmologyToolsPlugins)
set(PLUGIN_VERSION "1.0")

## Setup include directory paths
include_directories(
    ${PROJECT_SOURCE_DIR}/readers
    ${PROJECT_SOURCE_DIR}/filters
    ${PARAVIEW_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
    )

## Note we need to make sure that we put this include directory last so
## that the CosmologyTools halo-finder will be used and not the old out-dated
## version of the code that is available in ParaView.
include_directories(BEFORE ${PACKAGE_halofinder_INCLUDE_DIRS})

## Set the server-manager XML
set(SMXML
    resources/DensityProfile.xml
    resources/LANLHaloFinder.xml
    resources/vtkPVoronoiReader.xml
    )

## Set the server-manager sources
set(SMSRC
    filters/vtkCosmoDensityProfile.cxx
    filters/vtkPLANLHaloFinder.cxx
    readers/vtkPVoronoiReader.cxx
    )

## Add halo-finders to the project if that plug-in is built
set(RequiredLibs
    ${PACKAGE_halofinder_LIBRARIES}
    ${MPI_LIBRARIES}
    )

## Add the paraview plugin
ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION}
    GUI_RESOURCE_FILES resources/CosmologyTools.xml resources/vtkPVoronoiReaderGUI.xml
    SERVER_MANAGER_XML ${SMXML}
    SERVER_MANAGER_SOURCES ${SMSRC}
    )
target_link_libraries(${PLUGIN_NAME} ${RequiredLibs})
