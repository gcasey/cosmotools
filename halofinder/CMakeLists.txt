project(halofinder)

## Append this to the global include directories for this package
set(PACKAGE_halofinder_INCLUDE_DIRS
    ${PACKAGE_halofinder_INCLUDE_DIRS}
    ${PROJECT_SOURCE_DIR}
    CACHE INTERNAL "package include directories")
set(PACKAGE_halofinder_LIBRARIES
    ${PACKAGE_halofinder_LIBRARIES}
    fwdhalotracker
    CACHE INTERNAL "package libraries")

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Append to the include directories
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/cosmohalo)

## Add halofinder directories
add_subdirectory(cosmohalo)

## List of sources for the halo tracker
set(HALOTRACKER_SRC
    ForwardHaloTracker.cxx
    HaloDataInformation.cxx
    TemporalHaloInformation.cxx)

## Specify global source list for this project
#set(PACKAGE_halofinder_SOURCES
#    ${PACKAGE_halofinder_SOURCES}
#    ${HALOTRACKER_SRC}
#    CACHE INTERNAL "package sources")
GetPackageSources(halofinder ${PROJECT_SOURCE_DIR} "${HALOTRACKER_SRC}")

if(BUILD_SHARED_LIBS)
 add_library(fwdhalotracker SHARED ${HALOTRACKER_SRC})
else()
 add_library(fwdhalotracker STATIC ${HALOTRACKER_SRC})
endif()
target_link_libraries(fwdhalotracker
    ${PACKAGE_halofinder_LIBRARIES} ${MPI_C_LIBRARIES})
