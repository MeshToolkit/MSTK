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
pkg_check_modules(PC_PTScotch Quiet ptscotch)


# Search for include files

find_path(PTScotch_INCLUDE_DIR
  NAMES ptscotch.h
  HINTS ${PC_PTScotch_INCLUDE_DIRS}
  PATH_SUFFIXES parmetis)

if (NOT PTScotch_INCLUDE_DIR)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Cannot locate ptscotch.h")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIET)
      message(WARNING "Cannot locate ptscotch.h")
    endif ()
  endif ()
endif ()

set(PTScotch_INCLUDE_DIRS "${PTScotch_INCLUDE_DIR}")


# Search for libraries (scotch, ptscotch, ptscotcherr)

find_library(Scotch_LIBRARY
  NAMES scotch
  HINTS ${PC_PTScotch_LIBRARY_DIRS})

if (NOT Scotch_LIBRARY)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Can not locate Scotch library")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIET)
      message(WARNING "Cannot locate Scotch library")
    endif ()
  endif ()
endif ()

find_library(PTScotch_LIBRARY
  NAMES ptscotch
  HINTS ${PC_PTScotch_LIBRARY_DIRS})

if (NOT PTScotch_LIBRARY)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Can not locate PTScotch library")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIET)
      message(WARNING "Cannot locate PTScotch library")
    endif ()
  endif ()
endif ()

find_library(PTScotcherr_LIBRARY
  NAMES ptscotcherr
  HINTS ${PC_PTScotch_LIBRARY_DIRS})

if (NOT PTScotcherr_LIBRARY)
  if (PTScotch_FIND_REQUIRED)
    message(FATAL "Can not locate PTScotcherr library")
  else (PTScotch_FIND_REQUIRED)
    if (NOT PTScotch_FIND_QUIET)
      message(WARNING "Cannot locate PTScotcherr library")
    endif ()
  endif ()
endif ()

set(PTScotch_VERSION ${PC_PTScotch_VERSION})  # No guarantee

# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PTScotch
  FOUND_VAR PTScotch_FOUND
  REQUIRED_VARS
  Scotch_LIBRARY PTScotch_LIBRARY PTScotcherr_LIBRARY
  PTScotch_INCLUDE_DIR)


# Create PTScotch target and label scotch and ptscotcherr as dependencies

if (PTScotch_FOUND AND NOT TARGET MSTK::PTScotch)
  set(PTScotch_LIBRARIES MSTK::PTScotch)
  add_library(${PTScotch_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${PTScotch_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${PTScotch_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_PTScotch_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}")

  add_library(MSTK::Scotch UNKNOWN IMPORTED)
  set_target_properties(MSTK::Scotch PROPERTIES
    IMPORTED_LOCATION "${Scotch_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_PTScotch_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}")

  add_library(MSTK::PTScotcherr UNKNOWN IMPORTED)
  set_target_properties(MSTK::PTScotcherr PROPERTIES
    IMPORTED_LOCATION "${PTScotcherr_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_PTScotch_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${PTScotch_INCLUDE_DIR}")

  target_link_libraries(${PTScotch_LIBRARIES} INTERFACE MSTK::Scotch)
  target_link_libraries(${PTScotch_LIBRARIES} INTERFACE MSTK::PTScotcherr)
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

