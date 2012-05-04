project(interfaces)


## Add the interfaces binary directory to the list of include directories
## since that's where the mangled header used in the Fortran interface defined.
include_directories(${PROJECT_BINARY_DIR}/C)

message("Including: " ${PROJECT_BINARY_DIR}/C)

## We always build the C-interface
add_subdirectory(C)

option(BUILD_FORTRAN_INTERFACE "Build Fortran interface" OFF)
if( ${BUILD_FORTRAN_INTERFACE} )
  enable_language(Fortran)
  add_subdirectory(Fortran)
endif()

