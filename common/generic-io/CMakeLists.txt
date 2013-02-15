
## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
  get_filename_component(headerFile ${header} NAME)
  file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

## GenericIO source files
set(GENERIC_IO_SRC
    ${PROJECT_SOURCE_DIR}/generic-io/GenericIO.cxx
    )

## Append GenericIO sources to list of sources
set(COMMON_SRC
    ${COMMON_SRC}
    ${GENERIC_IO_SRC}
    CACHE INTERNAL ""
    )
GetPackageSources(common ${PROJECT_SOURCE_DIR}/generic-io "${SOURCES}")
