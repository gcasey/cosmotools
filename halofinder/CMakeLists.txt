project(halofinder)

option(BUILD_COSMO "Build the cosmo halo-finder." OFF)
if( ${BUILD_COSMO} )
  add_subdirectory(cosmohalo)
endif()
