# - Try to find the Qhull library in the system
# Once done this will define
#
# QHULL_FOUND -- boolean that indicates whether qhull was found
# QHULL_INCLUDE_DIR -- the include path for qhull
# QHULL_LIBRARIES -- the qhull libraries to link against

## Try to find the include directory
find_path(QHULL_INCLUDE_DIRS
 NAMES qhull_a.h
 PATHS /usr/include /usr/local/include /usr/include/qhull)

## Try to find the Qhull library
find_library(QHULL_LIBRARIES
 NAMES qhull
 PATHS /usr/lib64 /usr/lib)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    QHULL DEFAULT_MSG QHULL_INCLUDE_DIRS QHULL_LIBRARIES)
mark_as_advanced(QHULL_INCLUDE_DIRS)
mark_as_advanced(QHULL_LIBRARIES)
