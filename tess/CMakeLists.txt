project(tess)
message(STATUS "-- Building tess")

## Specify global include dirs and global libraries for this project
## such that, they can be used within other projects.
set(PACKAGE_tess_INCLUDE_DIRS
    ${PACKAGE_tess_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}/lib
    CACHE INTERNAL "package include directories")
set(PACKAGE_tess_LIBRARIES
    ${PACKAGE_tess_LIBRARIES}
    tess
    CACHE INTERNAL "package libraries")

## Move header files to the include directory in the build tree
file(GLOB HEADERS lib/*.h lib/*.hpp)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Find OpenGL library
find_package(OpenGL REQUIRED)
if( NOT OPENGL_FOUND )
 message(FATAL_ERROR "Could not find OpenGL!")
endif()
message(STATUS "-- Found OpenGL!")

## Find GLUT library
find_package(GLUT REQUIRED)
if( NOT GLUT_FOUND )
 message(FATAL_ERROR "Could not find Glut!")
endif()
message(STATUS "-- Glut is found!")

## Find DIY library
find_package(DIY REQUIRED)
if(NOT DIY_FOUND)
 message(FATAL_ERROR "Could not find DIY!")
endif()
message(STATUS "-- DIY is found!")

## Find Qhull
find_package(Qhull REQUIRED)
if(NOT QHULL_FOUND)
 message(FATAL_ERROR "Could not find Qhull!")
endif()
message(STATUS "-- Qhull is found!")

include_directories(
    ${DIY_INCLUDE_DIRS}
    ${QHULL_INCLUDE_DIRS}
    ${MPI_INCLUDE_DIRS}
    ${GLUT_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR}
    ${PROJECT_SOURCE_DIR}/lib)

## Required libraries for tess
set(RequiredLibs
    ${DIY_LIBRARIES}
    ${QHULL_LIBRARIES}
    ${MPI_LIBRARIES}
    ${GLUT_LIBRARIES}
    ${OPENGL_LIBRARIES}
    )

## TESS library sources
set(TESSLIB_SRC
  lib/ser_io.cpp
  lib/swap.cpp
  lib/tess.c
 )

## Set executables
set(TESSLIB_EXE
    lib/draw.cpp
    lib/dump.cpp
    driver/main.cpp
    )

## Build tess library
if(BUILD_SHARED_LIBS)
  add_library(tess SHARED ${TESSLIB_SRC})
else()
  add_library(tess STATIC ${TESSLIB_SRC})
endif()
target_link_libraries(tess ${RequiredLibs})

## Build executables
foreach(fileName ${TESSLIB_EXE})
  get_filename_component(myFile ${fileName} NAME)
  get_filename_component(mySuffix ${myFile} NAME_WE)
  add_executable(tess-${mySuffix} ${fileName})
  target_link_libraries(tess-${mySuffix} tess ${RequiredLibs})
endforeach()
