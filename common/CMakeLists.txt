project(commmon)

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## Configure CosmologyToolsDefinitions
configure_file(
    ${PROJECT_SOURCE_DIR}/CosmologyToolsDefinitions.h.in
    ${HEADER_INCLUDES_DIRECTORY}/CosmologyToolsDefinitions.h)
