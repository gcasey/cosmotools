cmake_minimum_required(VERSION 2.8)
project(pv-plugins)

## Find ParaView build and import cmake configurations, macros, etc.
find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})

## Set the plugin name and version
set(PLUGIN_NAME PVCosmologyToolsPlugins)
set(PLUGIN_VERSION "1.0")

include_directories(
    ${PROJECT_SOURCE_DIR}/readers
    ${PROJECT_SOURCE_DIR}/filters
    ${PARAVIEW_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
    ${PACKAGE_halofinder_INCLUDE_DIRS}
    )

## Set the server-manager XML
set(SMXML
    resources/DensityProfile.xml
    resources/LANLHaloFinder.xml
    )

## Set the server-manager sources
set(SMSRC
    filters/vtkCosmoDensityProfile.cxx
    filters/vtkPLANLHaloFinder.cxx
    )

## Add halo-finders to the project if that plug-in is built
set(RequiredLibs
    ${PACKAGE_halofinder_LIBRARIES})

## Append halo-finder XML definitions to the server-manager XML
#list(APPEND SMXML resources/LANLHaloFinder.xml)
## Append halo-finder sources to the server-manager XML
#list(APPEND SMSRC filters/vtkPLANLHaloFinder.cxx)


## Add the paraview plugin
ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION}
    GUI_RESOURCE_FILES resources/CosmologyTools.xml
    SERVER_MANAGER_XML ${SMXML}
    SERVER_MANAGER_SOURCES ${SMSRC}
    )
target_link_libraries(${PLUGIN_NAME} ${RequiredLibs})
