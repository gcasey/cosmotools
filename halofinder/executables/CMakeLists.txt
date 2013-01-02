include_directories(${PACKAGE_halofinder_INCLUDE_DIRS})

## list of executables
set(executables
    haloExchangeDIY.cxx)

## Setup required libraries
set(RequiredLibs
    ${MPI_C_LIBRARIES}
    ${DIY_LIBRARIES}
    ${PACKAGE_halofinder_LIBRARIES})

foreach(exeFile ${executables})
  get_filename_component(exe ${exeFile} NAME_WE)
  add_executable(${exe} ${exeFile})
  target_link_libraries(${exe} ${RequiredLibs})
endforeach()
