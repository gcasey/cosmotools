project(common)

## Specify global include dirs and global libraries for this project
## such that, they can be used within other projects.
set(PACKAGE_common_INCLUDE_DIRS
    ${PACKAGE_common_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "package include directories")
set(PACKAGE_common_LIBRARIES
    ${PACKAGE_common_LIBRARIES}
    common
    CACHE INTERNAL "package libraries")

set(SOURCES
    GenericIO.cxx
    )

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Specify global sources for this project
#set(PACKAGE_sfprobe_SOURCES
#    ${SOURCES}
#   CACHE INTERNAL "package sources")
GetPackageSources(common ${PROJECT_SOURCE_DIR} "${SOURCES}")

## Configure CosmologyToolsDefinitions
configure_file(
    ${PROJECT_SOURCE_DIR}/CosmologyToolsDefinitions.h.in
    ${HEADER_INCLUDES_DIRECTORY}/CosmologyToolsDefinitions.h)

if(BUILD_SHARED_LIBS)
 add_library(common SHARED ${SOURCES})
else()
 add_library(common STATIC ${SOURCES})
endif()
target_link_libraries(common ${MPI_C_LIBRARIES})

add_subdirectory(executables)
