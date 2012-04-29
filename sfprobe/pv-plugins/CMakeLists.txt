cmake_minimum_required(VERSION 2.8)
project(pv-plugins)

## Find ParaView build and import cmake configurations, macros, etc.
find_package(ParaView REQUIRED)
include(${PARAVIEW_USE_FILE})

## Set the plugin name and version
set(PLUGIN_NAME HACCTools)
set(PLUGIN_VERSION "1.0")

include_directories(
    ${PROJECT_SOURCE_DIR}/readers
    ${PROJECT_SOURCE_DIR}/filters
    ${PARAVIEW_INCLUDE_DIRS}
    ${VTK_INCLUDE_DIR}
    )

## Set the server-manager XML
set(SMXML
    resources/DensityProfile.xml
    )

## Set the server-manager sources
set(SMSRC
    filters/vtkCosmoDensityProfile.cxx
    )

## Add the paraview plugin
ADD_PARAVIEW_PLUGIN(${PLUGIN_NAME} ${PLUGIN_VERSION}
    GUI_RESOURCE_FILES resources/HACCTools.xml
    SERVER_MANAGER_XML ${SMXML}
    SERVER_MANAGER_SOURCES ${SMSRC}
    )
