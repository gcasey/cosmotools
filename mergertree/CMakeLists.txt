project(mergertree)

## Set all the merger-tree executables
set(MERGERTREE_EXE
    MergerTree.cpp)

foreach(fileName ${MERGERTREE_EXE})
  get_filename_component(prgm ${fileName} NAME_WE)
  add_executable(${prgm} ${fileName})
endforeach()


