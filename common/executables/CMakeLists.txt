include_directories(${PACKAGE_common_INCLUDE_DIRS})

## list of executables
set(executables
    GenericIOInspector.cxx
    )

## Setup required libraries
set(RequiredLibs
    ${MPI_C_LIBRARIES}
    ${PACKAGE_common_LIBRARIES}
    )

foreach(exeFile ${executables})
  get_filename_component(exe ${exeFile} NAME_WE)
  add_executable(${exe} ${exeFile})
  target_link_libraries(${exe} ${RequiredLibs})
endforeach()
