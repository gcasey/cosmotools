set(FSOURCES
    ModCosmologyTools.f90
    )

if(BUILD_SHARED_LIBS)
  add_library(fCosmologyToolsAPI SHARED ${FSOURCES})
else()
  add_library(fCosmologyToolsAPI STATIC ${FSOURCES})
endif()
target_link_libraries(fCosmologyToolsAPI cCosmologyToolsAPI ${MPI_LIBRARIES})
