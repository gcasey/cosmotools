set(FSOURCES
    ModCosmologyTools.f90
    )

## Add the MPI fortran include directory to the path,e.g., for mpi.mod
include_directories(${MPI_Fortran_INCLUDE_PATH})

## Set the list of required libraries
set(RequiredLibs
    cCosmologyToolsAPI
    ${MPI_C_LIBRARIES}
    ${MPI_Fortran_LIBRARIES}
    )

if(BUILD_SHARED_LIBS)
  add_library(fCosmologyToolsAPI SHARED ${FSOURCES})
else()
  add_library(fCosmologyToolsAPI STATIC ${FSOURCES})
endif()
target_link_libraries(fCosmologyToolsAPI ${RequiredLibs})
