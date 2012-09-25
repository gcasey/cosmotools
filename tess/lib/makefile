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

INCLUDE = -I$(HOME)/diy/include -I$(HOME)/software/qhull-2011.2/src/libqhull
LIBS = -L$(HOME)/diy/lib -ldiy -L$(HOME)/software/qhull-2011.2/lib -lqhullstatic -fopenmp

###########

ifeq ($(ARCH), MAC_OSX) # mac version

CCFLAGS = -g -Wall -Wextra
CCFLAGS += -DMAC_OSX $(TIMING) -fopenmp
GL_LIBS = -framework GLUT -framework OpenGL
CC=mpicc
CXX=g++

endif

###########

ifeq ($(ARCH), LINUX) # linux version

CCFLAGS = -g
CCFLAGS += -DLINUX $(TIMING)
GL_LIBS = -lglut -lGLU -lGL
CC=mpicc
CXX=g++

endif

###########

ifeq ($(ARCH), BGP) # BG/P version

CCFLAGS = -O3 -qarch=450d -qtune=450
#CCFLAGS = -g
CCFLAGS += -DBGP $(TIMING)
GL_LIBS =
CC=mpixlc_r
CXX=mpixlcxx_r
#CC=mpicxx
#CXX=mpicxx

endif

###########

TESS_OBJS = tess.o
DRAW_OBJS = draw.o ser_io.o swap.o
DUMP_OBJS = dump.o ser_io.o swap.o

.SUFFIXES: .c .cpp

.c.o:
	$(CC) -c $(CCFLAGS) $(INCLUDE) $<

.cpp.o:
	$(CXX) -c $(CCFLAGS) $(INCLUDE) $<

default: all

ifeq ($(ARCH), BGP)
all: libtess.a dump
else
all: libtess.a draw dump
endif

libtess.a: $(TESS_OBJS)
	rm -f $@
	ar -cru $@ $(TESS_OBJS)

draw: $(DRAW_OBJS) 
	$(CXX) -o draw $(DRAW_OBJS) $(GL_LIBS)

dump: $(DUMP_OBJS) 
	$(CXX) -o dump $(DUMP_OBJS)

clean:
	rm -f *.o 
