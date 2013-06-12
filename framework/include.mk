## Specify framework
FRAMEWORK_HOME := ${COSMOTOOLS_BASEDIR}/framework
FRAMEWORK_INCLUDES := ${FRAMEWORK_HOME}

## List framework sources
FRAMEWORK_SOURCES += ${FRAMEWORK_HOME}/InSituAlgorithm.cxx
FRAMEWORK_SOURCES += ${FRAMEWORK_HOME}/InSituAlgorithmInstantiator.cxx
FRAMEWORK_SOURCES += ${FRAMEWORK_HOME}/InSituAnalysisConfig.cxx
FRAMEWORK_SOURCES += ${FRAMEWORK_HOME}/InSituAnalysisManager.cxx
FRAMEWORK_SOURCES += ${FRAMEWORK_HOME}/LANLHaloFinderInSituAlgorithm.cxx
FRAMEWORK_SOURCES += ${FRAMEWORK_HOME}/SimulationParticles.cxx
 
## If DIY and QHull are defined, compile tess
ifdef ${USEDIY}
	ifdef ${USEQHULL}
		FRAMEWORK_SOURCES += ${FRAMEWORK_HOME}/TessInSituAlgorithm.cxx
	endif
endif