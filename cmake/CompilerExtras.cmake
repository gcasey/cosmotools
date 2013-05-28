set(CXXWARNINGS)
set(CWARNINGS)
if(CMAKE_COMPILER_IS_GNUCC)
  option(ENABLE_WARNINGS "Enable Warnings" OFF)
  if(${ENABLE_WARNINGS})
    set(CXXWARNINGS "-Wnon-virtual-dtor -Wno-long-long -ansi -Wcast-align -Wchar-subscripts -Wall -Wextra -Wpointer-arith -Wformat-security -Woverloaded-virtual -Wshadow -Wunused-parameter -fno-check-new -fno-common")
    set(CWARNINGS "-Wno-long-long -ansi -Wcast-align -Wchar-subscripts -Wall -Wextra -Wpointer-arith -Wformat-security -Wshadow -Wunused-parameter -fno-common")
  endif() # END if warnings are enabled
endif() # END if compiler is GNU

set(CMAKE_CXX_FLAGS ${CXXWARNINGS})
set(CMAKE_C_FLAGS ${CWARNINGS})
