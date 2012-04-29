project(sfprobe)

set(SOURCES
    ExtentPartitioner.cxx
    StructureFormationProbe.cxx
    ParallelStructureFormationProbe.cxx
    )

if(BUILD_SHARED_LIBS)
 add_library(sfprobe SHARED ${SOURCES})
else()
 add_library(sfprobe STATIC ${SOURCES})
endif()

add_executable(sfprobe.app sfprobe.cxx)
target_link_libraries(sfprobe.app sfprobe ${MPI_LIBRARIES})
