#-------------------------------------------------------------------------------
# Makefile  that illustrates how to link to CosmoTools 
#-------------------------------------------------------------------------------

## Setup path to DIY, used within COSMOTOOLS
DIYINCLUDEDIR=/opt/diy-revision-165/diy/include
DIYLIBDIR=/opt/diy-revision-165/diy/lib

## Setup path to COSMOTOOLS
COSMOTOOLSDIR=/home/gzagaris/Dev/CosmologyTools/Build/master
COSMOTOOLSINCLUDE=$(COSMOTOOLSDIR)/include
COSMOTOOLSLIBDIR=$(COSMOTOOLSDIR)/libs

## Setup includes and lib dependencies
INCLUDE=-I${DIYINCLUDEDIR} -I${COSMOTOOLSINCLUDE}

LIBS=-L${COSMOTOOLSLIBDIR} -lcosmotools -L${DIYLIBDIR} -ldiy

## Setup Compiler toolchain
CXX=mpicxx
CXXFLAGS= -g

hacc_particles_io.o: hacc_particles_io.h hacc_particles_io.cc
	$(CXX) $(INCLUDE) $(CXXFLAGS) -c hacc_particles_io.cc -o hacc_particles_io.o
    
libhaccio.a: hacc_particles_io.o
		ar rcs libhaccio.a hacc_particles_io.o
		ranlib libhaccio.a
    
ParticleReader: ParticleReader.cxx
	$(CXX) $(INCLUDE) $(CXXFLAGS) -o ParticleReader ParticleReader.cxx $(LIBS)
	
all: libhaccio.a ParticleReader

clean:
	rm -f *.o *.a ParticleReader