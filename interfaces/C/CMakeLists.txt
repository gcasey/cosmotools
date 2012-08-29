## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

set(SOURCES
 CosmologyTools.cxx
 )

include_directories(
    ${PACKAGE_framework_INCLUDE_DIRS}
    )

set(RequiredLibs
    ${PACKAGE_framework_LIBRARIES}
    ${MPI_LIBRARIES})

## Setup the fortran mangling
include(FortranCInterface)
FortranCInterface_HEADER(
    CosmologyToolsAPIMangling.h SYMBOLS
        cosmotools_fortran_initialize
        cosmotools_set_analysis_config
        cosmotools_set_domain_parameters
        cosmotools_set_particles
        cosmotools_coprocess
        cosmotools_finalize
        )

## Copy the mangled header to the main includes directory
file(COPY ${PROJECT_BINARY_DIR}/C/CosmologyToolsAPIMangling.h
      DESTINATION ${HEADER_INCLUDES_DIRECTORY} )

if(BUILD_SHARED_LIBS)
  add_library(cCosmologyToolsAPI SHARED ${SOURCES})
else()
  add_library(cCosmologyToolsAPI STATIC ${SOURCES})
endif()
target_link_libraries(cCosmologyToolsAPI ${RequiredLibs})
