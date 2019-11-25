# Copyright: 2019- Triad National Security, LLC
#
# PTScotch Find Module for MSTK
#
# PTScotch needs METIS; This module will try to find METIS as well and add it
# as a dependency
#
# Usage: To search a particular path you can specify the path in
# CMAKE_PREFIX_PATH, in the CMake variable PTScotch_DIR or environment
# variable PTScotch_ROOT
#
# Following variables are set:
# PTScotch_FOUND          (BOOL)   Flag indicating if PTScotch was found
# PTScotch_INCLUDE_DIRS   (PATH)   Path to PTScotch include files
# PTScotch_LIBRARY        (FILE)   PTScotch library (libzoltan.a, libzoltan.so)
# PTScotch_LIBRARIES      (LIST)   List of PTScotch targets
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_PTScotch QUIET ptscotch)


# Search for include files

find_path(PTScotch_INCLUDE_DIR
  NAMES ptscotch.h
  HINTS ${PC_PTScotch_INCLUDE_DIRS}
  PATHS ${PTScotch_DIR}
  PATH_SUFFIXES include)

if (NOT PTScotch_INCLUDE_DIR)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Cannot locate ptscotch.h")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIETLY)
      message(WARNING "Cannot locate ptscotch.h")
    endif ()
  endif ()
endif ()

set(PTScotch_INCLUDE_DIRS "${PTScotch_INCLUDE_DIR}")


# Search for libraries (scotch, ptscotch, ptscotcherr)

find_library(Scotch_LIBRARY
  NAMES scotch
  HINTS ${PC_PTScotch_LIBRARY_DIRS}
  PATHS ${PTScotch_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT Scotch_LIBRARY)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Can not locate Scotch library")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIETLY)
      message(WARNING "Cannot locate Scotch library")
    endif ()
  endif ()
endif ()

find_library(PTScotch_LIBRARY
  NAMES ptscotch
  HINTS ${PC_PTScotch_LIBRARY_DIRS}
  PATHS ${PTScotch_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT PTScotch_LIBRARY)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Can not locate PTScotch library")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIETLY)
      message(WARNING "Cannot locate PTScotch library")
    endif ()
  endif ()
endif ()

find_library(PTScotcherr_LIBRARY
  NAMES ptscotcherr
  HINTS ${PC_PTScotch_LIBRARY_DIRS}
  PATHS ${PTScotch_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT PTScotcherr_LIBRARY)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Can not locate PTScotcherr library")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIETLY)
      message(WARNING "Cannot locate PTScotcherr library")
    endif ()
  endif ()
endif ()

set(PTScotch_VERSION ${PC_PTScotch_VERSION})  # No guarantee

# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTScotch
  DEFAULT_MSG
  Scotch_LIBRARY PTScotch_LIBRARY PTScotcherr_LIBRARY
  PTScotch_INCLUDE_DIR)


# Create PTScotch target and label scotch and ptscotcherr as dependencies

if (PTScotch_FOUND AND NOT TARGET PTScotch::PTScotch)
  add_library(PTScotch::PTScotch UNKNOWN IMPORTED)
  set_target_properties(${PTScotch_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${PTScotch_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_PTScotch_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}")

  add_library(PTScotch::Scotch UNKNOWN IMPORTED)
  set_target_properties(PTScotch::Scotch PROPERTIES
    IMPORTED_LOCATION "${Scotch_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_PTScotch_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}")

  add_library(PTScotch::PTScotcherr UNKNOWN IMPORTED)
  set_target_properties(PTScotch::PTScotcherr PROPERTIES
    IMPORTED_LOCATION "${PTScotcherr_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_PTScotch_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}")

  set(PTScotch_LIBRARIES)
  list(APPEND PTScotch_LIBRARIES PTScotch::PTScotch PTScotch::Scotch PTScotch::PTScotcherr)
endif ()


# Hide these variables from the cache
mark_as_advanced(
  PTScotch_INCLUDE_DIR
  PTScotch_INCLUDE_DIRS
  PTScotch_LIBRARY
  PTScotch_LIBRARIES
  Scotch_LIBRARY
  PTScotcherr_LIBRARY
  )

