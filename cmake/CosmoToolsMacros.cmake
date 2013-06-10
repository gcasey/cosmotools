# - Given a list of headers, copy the to includes directory in the build tree.
# copy_headers(hdrs)
macro(copy_headers hdrs)
    foreach(hdr ${hdrs})
      file(COPY ${hdr} DESTINATION ${HEADER_INCLUDES_DIRECTORY})
    endforeach()
endmacro(copy_headers)

# - Builds a library for the associated code
# cosmotools_library(lib srccode)
# - Give the sources files, a library is compiled with the specified name. This
# library may be a shared library, a static library or an object library
# according to globally defined user-supplied parameters.
macro(cosmotools_library lib srccode)
    if(${BUILD_SINGLE_LIBRARY})
        add_library(${lib}_obj OBJECT ${srccode})
    elseif(BUILD_SHARED_LIBS)
        add_library(${lib} SHARED ${srccode})
    else()
        add_library(${lib} STATIC ${srccode})
    endif()
endmacro(cosmotools_library)
