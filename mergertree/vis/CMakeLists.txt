cmake_minimum_required(VERSION 2.8)

project(mtreeviz)

#find_package(VTK REQUIRED)
#include(${VTK_USE_FILE})

# VTK is used for infovis components.
find_package(VTK COMPONENTS
   vtkCommonDataModel
#   vtkViewsInfovis
   vtkRenderingCore
   vtkRenderingFreeTypeOpenGL
   vtkRenderingVolumeOpenGL
   NO_MODULE)

set_property(DIRECTORY APPEND PROPERTY COMPILE_DEFINITIONS ${VTK_DEFINITIONS})
include_directories(SYSTEM ${VTK_INCLUDE_DIRS})

add_executable(mtreeviz mtreeviz.cxx)
target_link_libraries(mtreeviz
    vtkCommonDataModel
#    vtkViewsInfovis
    vtkRenderingCore
    vtkRenderingFreeTypeOpenGL
    vtkRenderingVolumeOpenGL
    vtkIOImage)
