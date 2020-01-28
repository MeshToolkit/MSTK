# Copyright: 2019- Triad National Security, LLC
#
# UnitTest++ Find Module for MSTK
#
# Usage: To search a particular path you can specify the path in
# CMAKE_PREFIX_PATH, in the CMake variable UnitTest++_DIR or environment
# variable UnitTest++_ROOT
#
# Following variables are set:
# UnitTest++_FOUND          (BOOL)   Flag indicating if UnitTest++ was found
# UnitTest++_INCLUDE_DIRS   (PATH)   Path to UnitTest++ include files
# UnitTest++_LIBRARY        (FILE)   UnitTest++ library (libzoltan.a, libzoltan.so)
# UnitTest++_LIBRARIES      (LIST)   List of UnitTest++ targets (MSTK::UnitTest++)
# UnitTest++_ROOT           (PATH)   Top level directory where UnitTest++ is installed
# UnitTest++_DIR            (PATH)   Top level directory where UnitTest++ is installed
#
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_UnitTest++ QUIET unittest++)


# Search for include files

find_path(UnitTest++_INCLUDE_DIR
  NAMES unittest++.h UnitTest++.h
  HINTS ${PC_UnitTest++_INCLUDE_DIRS} ${UnitTest++_DIR}
  PATH_SUFFIXES include include/UnitTest++ include/unittest++)

if (NOT UnitTest++_INCLUDE_DIR)
  if (UnitTest++_FIND_REQUIRED)
    message(FATAL "Cannot locate unittest++.h or UnitTest++.h")
  else (UnitTest++_FIND_REQUIRED)
    if (NOT UnitTest++_FIND_QUIETLY)
      message(WARNING "Cannot locate unittest++.h or UnitTest++.h")
    endif ()
  endif ()
endif ()

set(UnitTest++_INCLUDE_DIRS "${UnitTest++_INCLUDE_DIR}")


# Search for libraries

find_library(UnitTest++_LIBRARY
  NAMES unittest++ UnitTest++
  HINTS ${PC_UnitTest++_LIBRARY_DIRS} ${UnitTest++_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT UnitTest++_LIBRARY)
  if (UnitTest++_FIND_REQUIRED)
    message(FATAL "Can not locate UnitTest++ library")
  else (UnitTest++_FIND_REQUIRED)
    if (NOT UnitTest++_FIND_QUIETLY)
      message(WARNING "Cannot locate UnitTest++ library")
    endif ()
  endif ()
endif ()

set(UnitTest++_VERSION PC_UnitTest++_VERSION})  # No guarantee

if (NOT UnitTest++_ROOT)
  set(UnitTest++_DIR "${UnitTest++_INCLUDE_DIR}/.." CACHE PATH "Top level dir of UnitTest++ installation" FORCE)
  set(UnitTest++_ROOT "${UnitTest++_INCLUDE_DIR}/.." CACHE PATH "Top level dir of UnitTest++ installation" FORCE)
endif ()


# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(UnitTest++
  DEFAULT_MSG
  UnitTest++_LIBRARY UnitTest++_INCLUDE_DIR UnitTest++_ROOT)

# find_package_handle_standard_args ignores case and sets PACKAGE_FOUND
if (NOT UnitTest++_FOUND AND UNITTEST++_FOUND)
  set(UnitTest++_FOUND ${UNITTEST++_FOUND})
endif ()

# Create UnitTest++ target

if (UnitTest++_FOUND AND NOT TARGET UnitTest++::UnitTest++)
  set(UnitTest++_LIBRARIES UnitTest++::UnitTest++)
  add_library(${UnitTest++_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${UnitTest++_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${UnitTest++_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_UnitTest++_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${UnitTest++_INCLUDE_DIR}")
endif ()


# Hide these variables from the cache
mark_as_advanced(
  UnitTest++_INCLUDE_DIR
  UnitTest++_INCLUDE_DIRS
  UnitTest++_LIBRARY
  UnitTest++_LIBRARIES
  )
