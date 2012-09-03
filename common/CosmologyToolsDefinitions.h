#ifndef COSMOLOGYTOOLSDEFINITIONS_H_
#define COSMOLOGYTOOLSDEFINITIONS_H_

// C/C++ required includes
#include <stdint.h>
#include <map>
#include <string>

/* Explicitly set precision for position/velocity
 * Behavior is controller by the user via a CMAKE build option.
 */
#ifdef DOUBLE_PRECISION_POSVEL
  typedef double POSVEL_T;
#else
  typedef float POSVEL_T;
#endif

/* Explicitly set precision for potential
 * Behavior is controller by the user via a CMAKE build option.
 */
#ifdef DOUBLE_PRECISION_POTENTIAL
  typedef double POTENTIAL_T;
#else
  typedef float POTENTIAL_T;
#endif

#ifdef DOUBLE_PRECISION_GRID
  typedef double GRID_T;
#else
  typedef float GRID_T;
#endif

/* Explicitly set whether to use 64-bit or 32-bit integers for ID types.
 * Behavior is controller by the user via a CMAKE build option.
 */
#ifdef USE_64_BIT_IDS
  typedef int64_t ID_T;
#else
  typedef int32_t ID_T;
#endif

/* Explicitly set the type for status/mask arrays
 * Behavior is hard-coded in this file.
 */
typedef int32_t STATUS_T;
typedef uint16_t MASK_T;

// Generic integer/floating point types

/* Set whether default floating type precision is double or single */
#ifdef DEFAULT_DOUBLE_PRECISION
  typedef double REAL;
#else
  typedef float REAL;
#endif

/* Set whether to use 64 or 32 bit by default for integer types */
#ifdef DEFAULT_64BIT_INTS
  typedef int64_t INTEGER;
#else
  typedef int32_t INTEGER;
#endif


/*
 * Define dictionary as key,value pair of strings. Used to store analysis
 * tool parameters
 */
typedef std::map<std::string,std::string> Dictionary;

#endif /* COSMOLOGYTOOLSDEFINITIONS_H_ */
