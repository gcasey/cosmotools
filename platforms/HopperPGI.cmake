## Configuration on Hopper using PGI compilers (default)

## setup compilers
set(BASE_COMPILER_DIR /opt/cray/xt-asyncpe/5.12/bin)
set(CMAKE_C_COMPILER ${BASE_COMPILER_DIR}/cc)
set(CMAKE_CXX_COMPILER ${BASE_COMPILER_DIR}/CC)
set(CMAKE_FORTRAN_COMPILER ${BASE_COMPILER_DIR}/ftn)

set(MPI_C_COMPILER ${CMAKE_C_COMPILER})
set(MPI_CXX_COMPILER ${CMAKE_CXX_COMPILER})
set(MPI_Fortran_COMPILER ${CMAKE_FORTRAN_COMPILER})

## set behavior for finding programs
set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
