## Specify base directory for the algorithms directory
ALGORITHMS_HOME := ${COSMOTOOLS_BASEDIR}/algorithms

## Append halofinder sources
include ${ALGORITHMS_HOME}/halofinder/include.mk
ALGORITHMS_INCLUDES += ${HALOFINDER_INCLUDES}
ALGORITHMS_SOURCES  += ${HALOFINDER_SOURCES}

