include_directories(${PACKAGE_sfprobe_INCLUDE_DIRS})

## list of executables
set(executables)

## list of libraries required to compile the executables
set(RequiredLibs
    sfprobe
    ${MPI_C_LIBRARIES}
    )

## Add VTK executables
if(${BUILD_VTK_PROGRAMS})

  set(headers
      probe.h
      )
  set(executables
      ${executables}
      probe.cxx
      )

  if(${ENABLE_DAX})
    set(executables
        ${executables}
        tbbProbe.cxx
        )
    set(cuda_executables
        cudaProbe.cu
        )
  endif()

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

if(Dax_Cuda_FOUND)
  foreach(exeFile ${cuda_executables})
    get_filename_component(exe ${exeFile} NAME_WE)
    cuda_add_executable(${exe} ${exeFile})
    target_link_libraries(${exe} ${RequiredLibs})
  endforeach()
endif()
