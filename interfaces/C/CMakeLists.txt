## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

set(SOURCES
 CosmologyToolsManager.cxx
 CosmologyTools.cxx
 )

include_directories(${PACKAGE_halofinder_INCLUDE_DIRS})

set(RequiredLibs
    ${PACKAGE_halofinder_LIBRARIES}
    ${MPI_LIBRARIES})

## Setup the fortran mangling
include(FortranCInterface)
FortranCInterface_HEADER(
    CosmologyToolsAPIMangling.h SYMBOLS
        CosmologyInit
        CosmologyFinit
        CosmologyEnableVis
        CosmologyDisableVis
        CosmologySetMemoryLayout
        CosmologySetTrackerFrequency
        CosmologySetParticles
        CosmologyFinalize
        )

if(BUILD_SHARED_LIBS)
  add_library(cCosmologyToolsAPI SHARED ${SOURCES})
else()
  add_library(cCosmologyToolsAPI STATIC ${SOURCES})
endif()
target_link_libraries(cCosmologyToolsAPI ${RequiredLibs})
