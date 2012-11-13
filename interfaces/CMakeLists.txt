project(interfaces)


## Add the interfaces binary directory to the list of include directories
## since that's where the mangled header used in the Fortran interface is
## defined.
include_directories(${PROJECT_BINARY_DIR}/C)

option(BUILD_FORTRAN_INTERFACE "Build Fortran interface" OFF)
if( ${BUILD_FORTRAN_INTERFACE} )
  enable_language(Fortran)
  add_subdirectory(C)
  add_subdirectory(Fortran)
else()
  add_subdirectory(C)  
endif()

