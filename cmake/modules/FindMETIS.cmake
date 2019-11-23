# Copyright: 2019- Triad National Security, LLC
#
# METIS Find Module for MSTK
#
# Usage: To search a particular path you can specify the path in
# CMAKE_PREFIX_PATH, in the CMake variable METIS_DIR or environment
# variable METIS_ROOT
#
# Following variables are set:
# METIS_FOUND          (BOOL)   Flag indicating if METIS was found
# METIS_INCLUDE_DIRS   (PATH)   Path to METIS include files
# METIS_LIBRARY        (FILE)   METIS library (libzoltan.a, libzoltan.so)
# METIS_LIBRARIES      (LIST)   List of METIS targets (MSTK::METIS)
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_METIS Quiet zoltan)


# Search for include files

find_path(METIS_INCLUDE_DIR
  NAMES metis.h
  HINTS ${PC_METIS_INCLUDE_DIRS}
  PATH_SUFFIXES zoltan)

if (NOT METIS_INCLUDE_DIR)
  if (METIS_FIND_REQUIRED)
    message(FATAL "Cannot locate zoltan.h")
  else (METIS_FIND_REQUIRED)
    if (NOT METIS_FIND_QUIET)
      message(WARNING "Cannot locate zoltan.h")
    endif ()
  endif ()
endif ()

set(METIS_INCLUDE_DIRS "${METIS_INCLUDE_DIR}")


# Search for libraries

find_library(METIS_LIBRARY
  NAMES metis
  HINTS ${PC_METIS_LIBRARY_DIRS})

if (NOT METIS_LIBRARY)
  if (METIS_FIND_REQUIRED)
    message(FATAL "Can not locate METIS library")
  else (METIS_FIND_REQUIRED)
    if (NOT METIS_FIND_QUIET)
      message(WARNING "Cannot locate METIS library")
    endif ()
  endif ()
endif ()

set(METIS_VERSION PC_METIS_VERSION})  # No guarantee


# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(METIS
  FOUND_VAR METIS_FOUND
  REQUIRED_VARS METIS_LIBRARY METIS_INCLUDE_DIR)


# Create METIS target

if (METIS_FOUND AND NOT TARGET MSTK::METIS)
  set(METIS_LIBRARIES MSTK::METIS)
  add_library(${METIS_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${METIS_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${METIS_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_METIS_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${METIS_INCLUDE_DIR}")
endif ()


# Hide these variables from the cache
mark_as_advanced(
  METIS_INCLUDE_DIR
  METIS_INCLUDE_DIRS
  METIS_LIBRARY
  METIS_LIBRARIES
  )
