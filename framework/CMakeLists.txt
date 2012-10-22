project(framework)

option(ENABLE_FRAMEWORK_STATISTICS "Enable Framework Statistics" OFF)
if(ENABLE_FRAMEWORK_STATISTICS)
  add_definitions(-DENABLESTATS)
endif()

## Specify global include dirs and global libraries for this project
## such that, they can be used within other projects.
set(PACKAGE_framework_INCLUDE_DIRS
    ${PACKAGE_framework_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "package include directories")
set(PACKAGE_framework_LIBRARIES
    ${PACKAGE_framework_LIBRARIES}
    framework
    CACHE INTERNAL "package libraries")

## Move header files to the include directory in the build tree
file(GLOB HEADERS *.h *.hpp)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Append to the include directories
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PACKAGE_halofinder_INCLUDE_DIRS}
    ${PACKAGE_tess_INCLUDE_DIRS}
    ${PACKAGE_tess_INCLUDE_DIRS}
    ${MPI_C_INCLUDE_PATH}
    )

## Set the required libraries for the framework
set(RequiredLibs
    ${PACKAGE_halofinder_LIBRARIES}
    ${PACKAGE_tess_LIBRARIES}
    ${PACKAGE_tess_EXTERNAL_LIBRARIES}
    ${MPI_C_LIBRARIES}
    )

## Set Framework code sources
set(FRAMEWORK_SRC
    AnalysisTool.cxx
    AnalysisToolInstantiator.cxx
    CosmologyToolsConfiguration.cxx
    CosmologyToolsManager.cxx
    LANLHaloFinderAnalysisTool.cxx
    SimulationParticles.cxx
    TessVoidFinderAnalysisTool.cxx
    )

## Specify global source list for this project
#set(PACKAGE_framework_SOURCES
#    ${FRAMEWORK_SRC}
#    CACHE INTERNAL "package sources")
GetPackageSources(framework ${PROJECT_SOURCE_DIR} "${FRAMEWORK_SRC}")

if(BUILD_SHARED_LIBS)
  add_library(framework SHARED ${FRAMEWORK_SRC})
else()
  add_library(framework STATIC ${FRAMEWORK_SRC})
endif()
target_link_libraries(framework
    ${PACKAGE_framework_LIBRARIES} ${RequiredLibs})

if(BUILD_TESTING)
 add_subdirectory(testing)
endif()
