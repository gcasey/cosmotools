## Specify base directory for the algorithms directory
ALGORITHMS_HOME := ${COSMOTOOLS_BASEDIR}/algorithms

## Append halofinder sources
include ${ALGORITHMS_HOME}/halofinder/include.mk
ALGORITHMS_INCLUDES += ${HALOFINDER_INCLUDES}
ALGORITHMS_SOURCES  += ${HALOFINDER_SOURCES}

## If DIY and QHull are defined, compile tess
ifdef ${USEDIY}
	ifdef ${USEQHULL}
		include ${ALGORITHMS_HOME}/tess/include.mk
		ALGORITHMS_INCLUDES += ${TESS_INCLUDES}
		ALGORITHMS_SOURCES	+= ${TESS_SOURCES}
	endif
endif