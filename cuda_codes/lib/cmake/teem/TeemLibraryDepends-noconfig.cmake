#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "teem" for configuration ""
set_property(TARGET teem APPEND PROPERTY IMPORTED_CONFIGURATIONS NOCONFIG)
set_target_properties(teem PROPERTIES
  IMPORTED_LOCATION_NOCONFIG "${_IMPORT_PREFIX}/lib/libteem.so.1.12.0"
  IMPORTED_SONAME_NOCONFIG "libteem.so.1"
  )

list(APPEND _IMPORT_CHECK_TARGETS teem )
list(APPEND _IMPORT_CHECK_FILES_FOR_teem "${_IMPORT_PREFIX}/lib/libteem.so.1.12.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
