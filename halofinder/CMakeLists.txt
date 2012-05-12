project(halofinder)

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

include_directories(${PROJECT_SOURCE_DIR})
option(BUILD_COSMO "Build the cosmo halo-finder." OFF)
if( ${BUILD_COSMO} )
  add_subdirectory(cosmohalo)
endif()
