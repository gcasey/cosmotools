## Set precision for particle position and velocity vectors
option(TYPE_POSVEL_DOUBLE
 "Use double-precision for particle position and velocity arrays" OFF)
if(${TYPE_POSVEL_DOUBLE})
 add_definitions(-DDOUBLE_PRECISION_POSVEL)
endif()

## Set whether to use double precision for potential
option(TYPE_POTENTIAL_DOUBLE
  "Use double-precision for particle potetnial arrays" OFF)
if(${TYPE_POTENTIAL_DOUBLE})
 add_definitions(-DDOUBLE_PRECISION_POTENTIAL)
endif()

option(TYPE_GRID_DOUBLE
    "Use double-precision for grid type" OFF)
if(${TYPE_GRID_DOUBLE})
 add_definitions(-DDOUBLE_PRECISION_GRID)
endif()

## Set whether to use 64-bit integer for IDs
option(TYPE_IDS_64BITS "Use 64 bit IDs for integers" OFF)
if(${TYPE_IDS_64BITS})
 add_definitions(-DUSE_64_BIT_IDS)
endif()

## Set whether to use double precision by default
option(TYPE_REAL_DOUBLE "Use double precision by default" OFF)
if(${TYPE_REAL_DOUBLE})
 add_definitions(-DDEFAULT_DOUBLE_PRECISION)
endif()

## Set whether to use 64-bit integers for all int types
option(TYPE_INT_64BITS "Use 64-bit integers by default" OFF)
if(${TYPE_INT_64BITS})
 add_definitions(-DDEFAULT_64BIT_INTS)
endif()
