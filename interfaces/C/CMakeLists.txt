set(SOURCES
 CosmologyToolsManager.cxx
 CosmologyTools.cxx
 )

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
        CosmologyFinalize
        )

if(BUILD_SHARED_LIBS)
  add_library(cCosmologyToolsAPI SHARED ${SOURCES})
else()
  add_library(cCosmologyToolsAPI STATIC ${SOURCES})
endif()
target_link_libraries(cCosmologyToolsAPI ${MPI_LIBRARIES})
