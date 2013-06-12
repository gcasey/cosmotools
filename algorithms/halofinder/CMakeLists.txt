project(halofinder)

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
copy_headers("${HEADERS}")


## Set cosmo sources
set(HALOFINDER_SRC
  BHTree.cxx
  ChainingMesh.cxx
  CosmoHaloFinder.cxx
  CosmoHaloFinderP.cxx
  FOFHaloProperties.cxx
  GridExchange.cxx
  HaloCenterFinder.cxx
  HaloFinderInput.cxx
  InitialExchange.cxx
  Message.cxx
  ParticleDistribute.cxx
  ParticleExchange.cxx
  Partition.cxx
  SODHalo.cxx
  SubHaloFinder.cxx
  Timer.cxx
  Timings.cxx
  bigchunk.cxx
  dims.cxx
  )

## Create the libraries
cosmotools_library(halofinder "${HALOFINDER_SRC}")
