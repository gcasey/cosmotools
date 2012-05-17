#----------------------------------------------------------------------------
#
# Tom Peterka
# Argonne National Laboratory
# 9700 S. Cass Ave.
# Argonne, IL 60439
# tpeterka@mcs.anl.gov
#
# (C) 2011 by Argonne National Laboratory.
# See COPYRIGHT in top-level directory.
#
#----------------------------------------------------------------------------

include ../user_defs.mk

INCLUDE = -I$(HOME)/diy/include -I../lib \
	-I$(HOME)/software/qhull-2011.2/src/libqhull

LIBS = 	-L../lib -ltess -L$(HOME)/diy/lib -ldiy \
	-L$(HOME)/software/qhull-2011.2/lib -lqhullstatic

###########

ifeq ($(ARCH), MAC_OSX) # mac version

CCFLAGS = -g
CCFLAGS += -DMAC_OSX
CXX=mpicxx

endif

###########

ifeq ($(ARCH), LINUX) # linux version

CCFLAGS = -g
CCFLAGS += -DLINUX
CXX=mpicxx

endif

###########

ifeq ($(ARCH), BGP) # BG/P version

CCFLAGS = -O3 -qarch=450d -qtune=450
#CCFLAGS = -g
CCFLAGS += -DBGP
CXX=mpixlcxx_r
#CXX=mpicxx

endif

###########

OBJS = main.o


.SUFFIXES: .cpp

.cpp.o:
	$(CXX) -c $(CCFLAGS) $(INCLUDE) $<

default: all

all: voronoi

voronoi: $(OBJS) ../lib/libtess.a
	$(CXX) -o voronoi $(OBJS) $(LIBS)
#	$(CXX) -o voronoi $(OBJS) $(LIBS) /home/tpeterka/software/valgrindobj/lib/valgrind/ppc32-bgp/valgrind_memcheck.o -Wl,-e,_start_valgrind
	cp voronoi ..

clean:
	rm -f *.o 
