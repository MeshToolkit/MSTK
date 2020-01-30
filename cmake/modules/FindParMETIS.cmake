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
# ParMETIS_LIBRARY        (FILE)   ParMETIS library (libparmetis.a, libparmetis.so)
# ParMETIS_LIBRARIES      (LIST)   List of ParMETIS targets (ParMETIS::ParMETIS)
# ParMETIS_ROOT           (PATH)   Top level directory where ParMETIS is installed
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_ParMETIS QUIET parmetis)


# Search for include files

find_path(ParMETIS_INCLUDE_DIR
  NAMES parmetis.h
  HINTS ${PC_ParMETIS_INCLUDE_DIRS}
  PATHS ${ParMETIS_DIR}
  PATH_SUFFIXES include)

if (NOT ParMETIS_INCLUDE_DIR)
  if (ParMETIS_FIND_REQUIRED)
    message(FATAL "Cannot locate parmetis.h")
  else ()
    if (NOT ParMETIS_FIND_QUIETLY)
      message(WARNING "Cannot locate parmetis.h")
    endif ()
  endif ()
endif ()

set(ParMETIS_INCLUDE_DIRS "${ParMETIS_INCLUDE_DIR}")


# Search for libraries

find_library(ParMETIS_LIBRARY
  NAMES parmetis
  HINTS ${PC_ParMETIS_LIBRARY_DIRS}
  PATHS ${ParMETIS_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT ParMETIS_LIBRARY)
  if (ParMETIS_FIND_REQUIRED)
    message(FATAL "Can not locate ParMETIS library")
  else (ParMETIS_FIND_REQUIRED)
    if (NOT ParMETIS_FIND_QUIETLY)
      message(WARNING "Cannot locate ParMETIS library")
    endif ()
  endif ()
endif ()

set(ParMETIS_VERSION ${PC_ParMETIS_VERSION})  # No guarantee

# Not sure if this is the right way to do it, but this is to help
# other upstream packages that attempt to find the ParMETIS package
# due to transitive dependencies
if (NOT ParMETIS_ROOT)
  get_filename_component(ParMETIS_ROOT "${ParMETIS_INCLUDE_DIR}/.." ABSOLUTE CACHE "Top level dir of ParMETIS installation" FORCE)
endif ()
if (NOT ParMETIS_DIR)
  get_filename_component(ParMETIS_DIR "${ParMETIS_INCLUDE_DIR}/.." ABSOLUTE CACHE "Top level dir of ParMETIS installation" FORCE)
endif ()



# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ParMETIS
  DEFAULT_MSG
  ParMETIS_LIBRARY ParMETIS_INCLUDE_DIR ParMETIS_ROOT)

# find_package_handle_standard_args ignores case and sets PACKAGE_FOUND
if (NOT ParMETIS_FOUND AND PARMETIS_FOUND)
  set(ParMETIS_FOUND ${PARMETIS_FOUND})
endif ()

# Create ParMETIS target

if (ParMETIS_FOUND AND NOT TARGET ParMETIS::ParMETIS)
  set(ParMETIS_LIBRARIES ParMETIS::ParMETIS)

  add_library(${ParMETIS_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${ParMETIS_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${ParMETIS_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_ParMETIS_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${ParMETIS_INCLUDE_DIR}")

  # ParMETIS depends on METIS. Attempt to find it
  
  # METIS does not (and likely will not) have a CMake config file but
  # one can dream
  
  set(metis_dir_save ${METIS_DIR})
  
  find_package(METIS QUIET CONFIG)
  
  if (NOT METIS_FOUND)
    set(METIS_DIR ${metis_dir_save})
    find_package(METIS QUIET REQUIRED MODULE)
  endif ()

  message(STATUS "ParMETIS_LIBRARIES ----> ${ParMETIS_LIBRARIES}")
  message(STATUS "METIS_LIBRARIES ----> ${METIS_LIBRARIES}")
  # Add METIS::METIS as a dependency of ParMETIS
  target_link_libraries(${ParMETIS_LIBRARIES} INTERFACE ${METIS_LIBRARIES})
endif ()
  

# Hide these variables from the cache
mark_as_advanced(
  ParMETIS_INCLUDE_DIR
  ParMETIS_INCLUDE_DIRS
  ParMETIS_LIBRARY
  ParMETIS_LIBRARIES
  )
