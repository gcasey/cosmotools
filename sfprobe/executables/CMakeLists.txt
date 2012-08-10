include_directories(${PACKAGE_sfprobe_INCLUDE_DIRS})

## list of executables
set(executables)

## list of libraries required to compile the executables
set(RequiredLibs
    sfprobe
    ${MPI_LIBRARIES}
    )

## Add VTK executables
if(${BUILD_VTK_PROGRAMS})
 set(executables
     ${executables}
     probe.cxx
     )

 set(RequiredLibs
     ${RequiredLibs}
     vtkCommonDataModel
     vtkFiltersCosmo
     vtkIOXML
     )
endif()

foreach(exeFile ${executables})
  get_filename_component(exe ${exeFile} NAME_WE)
  add_executable(${exe} ${exeFile})
  target_link_libraries(${exe} ${RequiredLibs})
endforeach()

