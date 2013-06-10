## MPI is a required dependency
find_package(MPI REQUIRED)
include_directories(
    ${MPI_C_INCLUDE_PATH})
add_definitions(-DOMPI_SKIP_MPICXX -DMPICH_SKIP_MPICXX)

## Enable DIY
option(ENABLE_DIY "Enable DIY" OFF)
if(${ENABLE_DIY})
  add_definitions(-DUSEDIY)
  find_package(DIY REQUIRED)
  include_directories(${DIY_INCLUDE_DIRS})
endif()

## Enable Qhull
option(ENABLE_QHULL "Enable Qhull" OFF)
if(${ENABLE_QHULL})
  add_definitions(-DUSEQHULL)
  find_package(Qhull REQUIRED)
  include_directories(${QHULL_INCLUDE_DIRS})
endif()

## Enable Dax
option(ENABLE_DAX "Enable DAX Toolkit" OFF)
if(${ENABLE_DAX})
  add_definitions(-DUSEDAX)
  find_package(Dax REQUIRED)
  DaxConfigureTBB(REQUIRED)
  DaxConfigureCuda()
endif()
