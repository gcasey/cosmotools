#!/bin/bash

nproc=8

bold=`tput bold`
normal=`tput sgr0`
baseDir=`pwd`
export CC=mpicc
export CXX=mpicxx


#==============================================================================
# STEP 1: Download & build DIY
echo "${bold}=== Download & Build DIY ===${normal}"
svn co -r178 https://svn.mcs.anl.gov/repos/diy/trunk diy-source

diySourceDir=$baseDir"/diy-source"
diyBuildDir=$baseDir"/diy-build"
mkdir $diyBuildDir

cd $diySourceDir
export MPICC=mpicc
export MPICXX=mpicxx
./configure --enable-fpic --disable-openmp --prefix=$diyBuildDir
make -j$nproc
make install

cd $baseDir

#==============================================================================
# STEP 2: Download & build Qhull
echo "${bold}=== Download & Build Qhull ===${normal}"

git clone https://github.com/gzagaris/gxzagas-qhull.git
mv gxzagas-qhull qhull-source

qhullSourceDir=$baseDir"/qhull-source"
qhullBuildDir=$baseDir"/qhull-build"
mkdir $qhullBuildDir

cd $qhullBuildDir
cmake -DCMAKE_INSTALL_PREFIX=$qhullBuildDir $qhullSourceDir
make -j$nproc
make install

cd $baseDir
