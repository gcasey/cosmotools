## Set list of tests for the framework
set(FRAMEWORK_TEST_SRC
  TestConfigurationParser.cxx
  )

## Set list of include directories
include_directories(
    ${PROJECT_SOURCE_DIR}
    ${PROJECT_SOURCE_DIR}/testing
    )

## Set list of required libraries
set(RequiredLibs
    ${PACKAGE_framework_LIBRARIES}
    ${MPI_C_LIBRARIES}
    )

## Compile all test sources
foreach(t ${FRAMEWORK_TEST_SRC})
  get_filename_component(myTest ${t} NAME_WE)
  add_executable( ${myTest} ${t} )
  target_link_libraries(${myTest} ${RequiredLibs} )
endforeach()

