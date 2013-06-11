##-----------------------------------------------------------------------------
# Top Level Makefile for CosmoTools
#
#------------------------------------------------------------------------------

## Get the basedirectory for CosmoTools
COSMOTOOLS_BASEDIR := $(shell pwd)

## Check if an object directory is defined in the environment
## in which case cosmotools will be compiled in the specified 
## directory.
ifndef ${COSMOTOOLS_OBJDIR}
	COSMOTOOLS_OBJDIR = ${COSMOTOOLS_BASEDIR}/cosmotools.build
endif

## If an mpicxx compiler is prescribed in the environment, use that.
## Otherwise, default to mpicxx being the user's path.
ifndef ${COSMOTOOLS_MPICXX}
	COSMOTOOLS_MPICXX := mpicxx
endif

## Propagate compile flags from the environment
ifndef ${COSMOTOOLS_CXXFLAGS}
	COSMOTOOLS_CXXFLAGS := -O3
endif

## include subdirectories
include common/include.mk
include algorithms/include.mk
#include framework/include.mk

## Setup the cosmotools includes
COSMOTOOLS_INCLUDES += ${COMMON_INCLUDES} ${ALGORITHMS_INCLUDES}

## Setup the cosmotools sources
COSMOTOOLS_SOURCES += ${COMMON_SOURCES} ${ALGORITHMS_SOURCES}

## Get list of object targets
OBJECTS = $(COSMOTOOLS_SOURCES:.cxx=.o)

##-----------------------------------------------------------------
##								T A	R	G	E	T	S
##-----------------------------------------------------------------

default: all

all: ${COSMOTOOLS_OBJDIR}/libcosmotools.a

## Target to make build directory
$(COSMOTOOLS_OBJDIR):
	mkdir -p $(COSMOTOOLS_OBJDIR)
	
%.o: %.cxx | $(COSMOTOOLS_OBJDIR)
	${COSMOTOOLS_MPICXX} ${COSMOTOOLS_INCLUDES} ${COSMOTOOLS_CXXFLAGS} -c -o $@ $<
	
$(COSMOTOOLS_OBJDIR)/libcosmotools.a: $(COSMOTOOLS_OBJDIR)/libcosmotools.a($(OBJECTS))
	ranlib $@
	
clean:
	-rm -rf $(COSMOTOOLS_OBJDIR)