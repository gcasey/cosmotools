##-----------------------------------------------------------------------------
# Top level Makefile for CosmoTools
#
#------------------------------------------------------------------------------

## Get the basedirectory for CosmoTools
COSMOTOOLS_BASEDIR := $(shell pwd)

# Used to get bold text in echo statements
BOLD		= "\033[1m"

# End bold text
NBOLD		= "\033[0m"

# Use this macro to print info to the console
ECHO		= /bin/echo -e

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
	COSMOTOOLS_CXXFLAGS := -O2 -Wall
endif

## include subdirectories
include common/include.mk
include algorithms/include.mk
#include framework/include.mk

## Setup the cosmotools includes
COSMOTOOLS_INCLUDES += -I${GENERIC_IO_INCLUDES}
COSMOTOOLS_INCLUDES += -I${COMMON_INCLUDES} 
COSMOTOOLS_INCLUDES += -I${ALGORITHMS_INCLUDES}

## Setup the cosmotools sources
COSMOTOOLS_SOURCES += ${COMMON_SOURCES} 
COSMOTOOLS_SOURCES += ${ALGORITHMS_SOURCES}

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
	$(ECHO) -n $(BOLD) " ** " $< " ** \n" $(NBOLD) 
	${COSMOTOOLS_MPICXX} ${COSMOTOOLS_INCLUDES} ${COSMOTOOLS_CXXFLAGS} -c $< -o $@
	
$(COSMOTOOLS_OBJDIR)/libcosmotools.a: $(COSMOTOOLS_OBJDIR)/libcosmotools.a($(OBJECTS))
	$(ECHO) -n $(BOLD) " ** " $@ " ** \n" $(NBOLD) 
	ranlib $@
	
clean:
	-rm -rf $(COSMOTOOLS_OBJDIR)