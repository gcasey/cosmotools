project(common)

set(COMMON_SRC
    MPIUtilities.cxx
    TaskTimer.cxx
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


if(BUILD_SHARED_LIBS)
 add_library(common SHARED ${COMMON_SRC})
else()
 add_library(common STATIC ${COMMON_SRC})
endif()
target_link_libraries(common ${MPI_C_LIBRARIES})

add_subdirectory(executables)
