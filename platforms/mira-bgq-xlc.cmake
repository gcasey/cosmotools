## Configuration on Hopper using PGI compilers (default)

## setup compilers
set(BASE_COMPILER_DIR /bgsys/drivers/ppcfloor/comm/xl/bin)
set(CMAKE_C_COMPILER ${BASE_COMPILER_DIR}/mpixlc_r)
set(CMAKE_CXX_COMPILER ${BASE_COMPILER_DIR}/mpixlcxx_r)
set(CMAKE_FORTRAN_COMPILER ${BASE_COMPILER_DIR}/mpixlf90_r)

set(GNUFLAGS "-O3")
set(CMAKE_CXX_FLAGS "${GNUFLAGS}")
set(CMAKE_C_FLAGS "${GNUFLAGS}")

set(CMAKE_FIND_ROOT_PATH
    	/bgsys/drivers/ppcfloor/comm/xl/)

set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
