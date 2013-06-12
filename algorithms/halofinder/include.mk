## Specify halofinder
HALOFINDER_HOME := ${ALGORITHMS_HOME}/halofinder
HALOFINDER_INCLUDES := ${HALOFINDER_HOME}

## List of halofinder sources
HALOFINDER_SOURCES := $(wildcard $(HALOFINDER_HOME)/*.cxx)