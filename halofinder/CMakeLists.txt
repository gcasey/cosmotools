project(halofinder)

## Append this to the global include directories for this package
set(PACKAGE_halofinder_INCLUDE_DIRS
    ${PACKAGE_halofinder_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "package include directories")

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Append to the include directories
include_directories(${PROJECT_SOURCE_DIR})

## Add halofinder directories
add_subdirectory(cosmohalo)
