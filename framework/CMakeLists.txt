project(framework)

## Move header files to the include directory in the build tree
file(GLOB HEADERS *.h *.hpp)
copy_headers("${HEADERS}")


## Set Framework code sources
set(FRAMEWORK_SRC
    InSituAlgorithm.cxx
    InSituAlgorithmInstantiator.cxx
    InSituAnalysisConfig.cxx
    InSituAnalysisManager.cxx
    LANLHaloFinderInSituAlgorithm.cxx
    SimulationParticles.cxx
    )

if(ENABLE_DIY AND ENABLE_QHULL)
  list(APPEND FRAMEWORK_SRC TessInSituAlgorith.cxx)
endif()


## Create the framework library
cosmotools_library(framework "${FRAMEWORK_SRC}")
