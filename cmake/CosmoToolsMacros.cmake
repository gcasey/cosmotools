# - Given a list of sources with relative path, get the list with absolute path
# GetPackageSources(pkg basepath srclist)
# - macro which takes a source list and converts it to
#
# The following variables are expected to have been set:
#
#   PACKAGE_${pkg}_SOURCES



macro(GetPackageSources pkg basepath srclist)
#    message(STATUS "Setting package sources for ${pkg}")
    foreach(src ${srclist})
#      message(STATUS "setting source: ${basepath}/${src}")
      set(PACKAGE_${pkg}_SOURCES
        ${PACKAGE_${pkg}_SOURCES}
        ${basepath}/${src}
        CACHE INTERNAL "${pkg} sources")
    endforeach()
endmacro(GetPackageSources)
