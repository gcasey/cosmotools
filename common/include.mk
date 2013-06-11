## Specify the include path of this directory
COMMON_HOME := ${COSMOTOOLS_BASEDIR}/common
COMMON_INCLUDES := ${COMMON_HOME}

## List of sources to be compiled under this directory
COMMON_SOURCES += ${COMMON_HOME}/MPIUtilities.cxx
COMMON_SOURCES += ${COMMON_HOME}/TaskTimer.cxx