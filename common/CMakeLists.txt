project(common)

## Sources in common
set(COMMON_SRC
    MPIUtilities.cxx
    TaskTimer.cxx
    )

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
copy_headers("${HEADERS}")

## Configure CosmoToolsDefinitions and put it in the HEADER_INCLUDES_DIR
configure_file(
    ${PROJECT_SOURCE_DIR}/CosmoToolsDefinitions.h.in
    ${HEADER_INCLUDES_DIRECTORY}/CosmoToolsDefinitions.h)


cosmotools_library(common "${COMMON_SRC}")
