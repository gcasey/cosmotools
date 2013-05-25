project(sqlite)

## SQLITE source files
set(SQLITE_SRC
    shell.c
    sqlite3.c
    )

## Move header files to include directory in the build tree
file(GLOB HEADERS *.h)
foreach(header ${HEADERS})
    get_filename_component(headerFile ${header} NAME)
    file(COPY ${header} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
endforeach()

if(BUILD_SHARED_LIBS)
    add_library(sqlite SHARED ${SQLITE_SRC})
else()
    add_library(sqlite STATIC ${SQLITE_SRC})
endif()
