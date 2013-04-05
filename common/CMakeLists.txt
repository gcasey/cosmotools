project(common)

## Specify global include dirs and global libraries for this project
## such that, they can be used within other projects.
set(PACKAGE_common_INCLUDE_DIRS
    ${PACKAGE_common_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/generic-io
    CACHE INTERNAL "package include directories")
set(PACKAGE_common_LIBRARIES
    ${PACKAGE_common_LIBRARIES}
    common
    CACHE INTERNAL "package libraries")

set(COMMON_SRC
    MPIUtilities.cxx
    CACHE INTERNAL "package sources")

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Configure CosmologyToolsDefinitions
configure_file(
    ${PROJECT_SOURCE_DIR}/CosmologyToolsDefinitions.h.in
    ${HEADER_INCLUDES_DIRECTORY}/CosmologyToolsDefinitions.h)


## Add generic-io
add_subdirectory(generic-io)

GetPackageSources(common ${PROJECT_SOURCE_DIR} "${COMMON_SRC}")

if(BUILD_SHARED_LIBS)
 add_library(common SHARED ${COMMON_SRC})
else()
 add_library(common STATIC ${COMMON_SRC})
endif()
target_link_libraries(common ${MPI_C_LIBRARIES})

add_subdirectory(executables)
