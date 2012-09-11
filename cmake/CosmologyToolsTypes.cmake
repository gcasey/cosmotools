## Set precision for particle position and velocity vectors
option(TYPE_POSVEL_DOUBLE
 "Use double-precision for particle position and velocity arrays" OFF)
if(${TYPE_POSVEL_DOUBLE})
 add_definitions(-DTYPE_POSVEL_DOUBLE)
endif()

## Set whether to use double precision for potential
option(TYPE_POTENTIAL_DOUBLE
  "Use double-precision for particle potetnial arrays" OFF)
if(${TYPE_POTENTIAL_DOUBLE})
 add_definitions(-DTYPE_POTENTIAL_DOUBLE)
endif()

option(TYPE_GRID_DOUBLE
    "Use double-precision for grid type" OFF)
if(${TYPE_GRID_DOUBLE})
 add_definitions(-DTYPE_GRID_DOUBLE)
endif()

## Set whether to use 64-bit integer for IDs
option(TYPE_IDS_64BITS "Use 64 bit IDs for integers" OFF)
if(${TYPE_IDS_64BITS})
 add_definitions(-DTYPE_IDS_64BITS)
endif()

## Set whether to use double precision by default
option(TYPE_REAL_DOUBLE "Use double precision by default" OFF)
if(${TYPE_REAL_DOUBLE})
 add_definitions(-DTYPE_REAL_DOUBLE)
endif()

## Set whether to use 64-bit integers for all int types
option(TYPE_INT_64BITS "Use 64-bit integers by default" OFF)
if(${TYPE_INT_64BITS})
 add_definitions(-DTYPE_INT_64BITS)
endif()
