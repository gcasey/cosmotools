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

  set(shared_srcs
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
        cudaProbe.cxx
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
  add_executable(${exe} ${exeFile} ${shared_srcs})
  target_link_libraries(${exe} ${RequiredLibs})
endforeach()

if(Dax_Cuda_FOUND)
  foreach(exeFile ${cuda_executables})
    get_filename_component(exe ${exeFile} NAME_WE)
    add_executable(${exe} ${exeFile} ${shared_srcs})
    set_target_properties(${exe} PROPERTIES COMPILE_FLAGS -DDAX_CUDA)
    target_link_libraries(${exe} cuda_sfprobe ${MPI_C_LIBRARIES} ${RequiredLibs})
  endforeach()
endif()
