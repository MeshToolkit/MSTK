# Copyright: 2019- Triad National Security, LLC
#
# ParMETIS Find Module for MSTK
#
# ParMETIS needs METIS; This module will try to find METIS as well and add it
# as a dependency
#
# Usage: To search a particular path you can specify the path in
# CMAKE_PREFIX_PATH, in the CMake variable ParMETIS_DIR or environment
# variable ParMETIS_ROOT
#
# Following variables are set:
# ParMETIS_FOUND          (BOOL)   Flag indicating if ParMETIS was found
# ParMETIS_INCLUDE_DIRS   (PATH)   Path to ParMETIS include files
# ParMETIS_LIBRARY        (FILE)   ParMETIS library (libzoltan.a, libzoltan.so)
# ParMETIS_LIBRARIES      (LIST)   List of ParMETIS targets (MSTK::ParMETIS)
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_ParMETIS Quiet parmetis)


# Search for include files

find_path(ParMETIS_INCLUDE_DIR
  NAMES parmetis.h
  HINTS ${PC_ParMETIS_INCLUDE_DIRS}
  PATH_SUFFIXES parmetis)

if (NOT ParMETIS_INCLUDE_DIR)
  if (ParMETIS_FIND_REQUIRED)
    message(FATAL "Cannot locate parmetis.h")
  else (ParMETIS_FIND_REQUIRED)
    if (NOT ParMETIS_FIND_QUIET)
      message(WARNING "Cannot locate parmetis.h")
    endif ()
  endif ()
endif ()

set(ParMETIS_INCLUDE_DIRS "${ParMETIS_INCLUDE_DIR}")


# Search for libraries

find_library(ParMETIS_LIBRARY
  NAMES parmetis
  HINTS ${PC_ParMETIS_LIBRARY_DIRS})

if (NOT ParMETIS_LIBRARY)
  if (ParMETIS_FIND_REQUIRED)
    message(FATAL "Can not locate ParMETIS library")
  else (ParMETIS_FIND_REQUIRED)
    if (NOT ParMETIS_FIND_QUIET)
      message(WARNING "Cannot locate ParMETIS library")
    endif ()
  endif ()
endif ()

set(ParMETIS_VERSION ${PC_ParMETIS_VERSION})  # No guarantee


# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS
  FOUND_VAR ParMETIS_FOUND
  REQUIRED_VARS ParMETIS_LIBRARY ParMETIS_INCLUDE_DIR)


# Create ParMETIS target

if (ParMETIS_FOUND AND NOT TARGET MSTK::ParMETIS)
  set(ParMETIS_LIBRARIES MSTK::ParMETIS)
  add_library(${ParMETIS_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${ParMETIS_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${ParMETIS_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_ParMETIS_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${ParMETIS_INCLUDE_DIR}")
endif ()


# ParMETIS depends on METIS. Attempt to find it

# METIS does not (and likely will not) have a CMake config file but
# one can dream

find_package(METIS CONFIG)
if (NOT METIS_FOUND)
  find_package(METIS REQUIRED MODULE)
endif ()

# Add METIS::METIS as a dependency of ParMETIS
target_link_libraries(${ParMETIS_LIBRARIES} INTERFACE ${METIS_LIBRARIES})
target_include_directories(${ParMETIS_LIBRARIES} INTERFACE ${METIS_INCLUDE_DIRS})


# Hide these variables from the cache
mark_as_advanced(
  ParMETIS_INCLUDE_DIR
  ParMETIS_INCLUDE_DIRS
  ParMETIS_LIBRARY
  ParMETIS_LIBRARIES
  )
