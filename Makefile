##-----------------------------------------------------------------------------
# Top level Makefile for CosmoTools
#
#------------------------------------------------------------------------------

## Get the basedirectory for CosmoTools
COSMOTOOLS_BASEDIR := $(shell pwd)

## Make modules path
MAKE_MODULES_PATH := $(COSMOTOOLS_BASEDIR)/make

## Load the GNU Make Standard Library (GMSL)
include $(MAKE_MODULES_PATH)/gmsl

## Global counter
COUNTER=0

# Use this macro to print info to the console
ECHO := @/bin/echo 

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

## If there is a user-supplied EXTERNAL_ALGORITHMS directory, use that
ifdef ${EXTERNAL_ALGORITHMS_DIRECTORY}
  include ${EXTERNAL_ALGORITHMS_DIRECTORY}/include.mk
endif

include framework/include.mk

## Setup the cosmotools includes
COSMOTOOLS_INCLUDES += -I${GENERIC_IO_INCLUDES}
COSMOTOOLS_INCLUDES += -I${COMMON_INCLUDES} 
COSMOTOOLS_INCLUDES += -I${ALGORITHMS_INCLUDES}
COSMOTOOLS_INCLUDES += -I${FRAMEWORK_INCLUDES}

## Setup the cosmotools sources
COSMOTOOLS_SOURCES += ${COMMON_SOURCES} 
COSMOTOOLS_SOURCES += ${ALGORITHMS_SOURCES}
COSMOTOOLS_SOURCES += ${FRAMEWORK_SOURCES}

## Get list of object targets
OBJECTS := $(COSMOTOOLS_SOURCES:.cxx=.o)
NOBJECTS := $(call length,$(OBJECTS))

##-----------------------------------------------------------------
##				  TARGETS
##-----------------------------------------------------------------

default: all

all: ${COSMOTOOLS_OBJDIR}/libcosmotools.a

## Target to make build directory
$(COSMOTOOLS_OBJDIR):
	$(ECHO) "==============================================================="
	$(ECHO) "     Building CosmoTools...                                    "
	$(ECHO) "==============================================================="
	@mkdir -p $(COSMOTOOLS_OBJDIR)
	
%.o: %.cxx | $(COSMOTOOLS_OBJDIR)
	$(eval COUNTER := $(call plus,$(COUNTER),1))
	$(ECHO) "[" $(COUNTER) "/" $(NOBJECTS) "]" $<   
	@${COSMOTOOLS_MPICXX} ${COSMOTOOLS_INCLUDES} ${COSMOTOOLS_CXXFLAGS} -c $< -o $@
	
$(COSMOTOOLS_OBJDIR)/libcosmotools.a: $(COSMOTOOLS_OBJDIR)/libcosmotools.a($(OBJECTS))
	$(ECHO) " *** " $@   
	@ranlib $@
	$(ECHO) "==============================================================="
	$(ECHO) "     Finished CosmoTools build                                 "
	$(ECHO) "==============================================================="
	
clean:
	-rm -rf $(COSMOTOOLS_OBJDIR)
