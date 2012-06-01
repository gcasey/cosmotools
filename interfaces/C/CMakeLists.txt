## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

set(SOURCES
 CosmologyTools.cxx
 CosmologyToolsManager.cxx
 SimulationParticles.cxx
 )

include_directories(
    ${PACKAGE_halofinder_INCLUDE_DIRS})

set(RequiredLibs
    ${PACKAGE_halofinder_LIBRARIES})

## find paraview
if( ${ENABLE_COVIS} )
    list(APPEND SOURCES ParaViewInSituInterface.cxx)

    find_package(ParaView REQUIRED)
    include(${PARAVIEW_USE_FILE})

    include_directories(
        ${VTK_INCLUDE_DIRS}
        ${PARAVIEW_INCLUDE_DIRS}
        )

   list(APPEND RequiredLibs ${PARAVIEW_LIBRARIES})
   list(APPEND RequiredLibs ${VTK_LIBRARIES})
   list(APPEND RequiredLibs ${MPI_LIBRARIES})
endif()

## Setup the fortran mangling
include(FortranCInterface)
FortranCInterface_HEADER(
    CosmologyToolsAPIMangling.h SYMBOLS
        CosmologyDisableVis
        CosmologyEnableVis
        CosmologyFinalize
        CosmologyFindHalos
        CosmologyFinit
        CosmologyGetHaloIds
        CosmologyGetSubHaloIds
        CosmologyInit
        CosmologySetAnalysisTimeSteps
        CosmologySetHaloFinder
        CosmologySetMemoryLayout
        CosmologySetParticles
        CosmologySetTrackerFrequency
        CosmologySetVisualizationFrequency
        CosmologyTrackHalos
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
