# Copyright: 2019- Triad National Security, LLC
#
# Zoltan Find Module for MSTK
#
# Usage: To search a particular path you can specify the path in
# CMAKE_PREFIX_PATH, in the CMake variable Zoltan_DIR or environment
# variable Zoltan_ROOT
#
# Following variables are set:
# Zoltan_FOUND          (BOOL)   Flag indicating if Zoltan was found
# Zoltan_INCLUDE_DIRS   (PATH)   Path to Zoltan include files
# Zoltan_LIBRARY        (FILE)   Zoltan library (libzoltan.a, libzoltan.so)
# Zoltan_LIBRARIES      (LIST)   List of Zoltan targets (MSTK::Zoltan)
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_Zoltan Quiet zoltan)


# Search for include files

find_path(Zoltan_INCLUDE_DIR
  NAMES zoltan.h
  HINTS ${PC_Zoltan_INCLUDE_DIRS} ${Zoltan_DIR} $ENV{Zoltan_ROOT}
  PATH_SUFFIXES zoltan)

if (NOT Zoltan_INCLUDE_DIR)
  if (Zoltan_FIND_REQUIRED)
    message(FATAL "Cannot locate zoltan.h")
  else (Zoltan_FIND_REQUIRED)
    if (NOT Zoltan_FIND_QUIET)
      message(WARNING "Cannot locate zoltan.h")
    endif ()
  endif ()
endif ()

set(Zoltan_INCLUDE_DIRS "${Zoltan_INCLUDE_DIR}")


# Search for libraries

find_library(Zoltan_LIBRARY
  NAMES zoltan
  HINTS ${PC_Zoltan_LIBRARY_DIRS} ${Zoltan_DIR} $ENV{Zoltan_ROOT})

if (NOT Zoltan_LIBRARY)
  if (Zoltan_FIND_REQUIRED)
    message(FATAL "Can not locate Zoltan library")
  else (Zoltan_FIND_REQUIRED)
    if (NOT Zoltan_FIND_QUIET)
      message(WARNING "Cannot locate Zoltan library")
    endif ()
  endif ()
endif ()

set(Zoltan_VERSION PC_Zoltan_VERSION})  # No guarantee


# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Zoltan
  FOUND_VAR Zoltan_FOUND
  REQUIRED_VARS Zoltan_LIBRARY Zoltan_INCLUDE_DIR)


# Create Zoltan target

if (Zoltan_FOUND AND NOT TARGET MSTK::Zoltan)
  set(Zoltan_LIBRARIES MSTK::Zoltan)
  add_library(${Zoltan_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${Zoltan_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${Zoltan_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_Zoltan_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${Zoltan_INCLUDE_DIR}")
endif ()


# Don't know how to automatically determine if Zoltan needs ParMETIS
# or PTScotch So, these are controlled by the variables
# Zoltan_NEEDS_ParMETIS and Zoltan_NEEDS_PTSCOTCH

if (Zoltan_NEEDS_ParMETIS)
  # ParMETIS does not (and likely will not) have a CMake config file
  # but one can dream

  find_package(ParMETIS CONFIG)
  if (NOT ParMETIS_FOUND)
    find_package(ParMETIS REQUIRED MODULE)
  endif ()

  # Add MSTK::ParMETIS as a dependency of Zoltan
  target_link_libraries(${Zoltan_LIBRARIES} INTERFACE ${ParMETIS_LIBRARIES})
  target_include_directories(${Zoltan_LIBRARIES} INTERFACE ${ParMETIS_INCLUDE_DIRS})
endif (Zoltan_NEEDS_ParMETIS)


if (Zoltan_NEEDS_PTScotch)
  # PTScotch does not have a CMake config file but one can dream

  find_package(PTScotch CONFIG)
  if (NOT PTScotch_FOUND)
    find_package(PTScotch REQUIRED MODULE)
  endif ()

  # Add MSTK::PTScotch as a dependency of Zoltan
  target_link_libraries(${Zoltan_LIBRARIES} INTERFACE ${PTScotch_LIBRARIES})
  target_include_directories(${Zoltan_LIBRARIES} INTERFACE ${PTScotch_INCLUDE_DIRS})
endif (Zoltan_NEEDS_PTScotch)


# Hide these variables from the cache
mark_as_advanced(
  Zoltan_INCLUDE_DIR
  Zoltan_INCLUDE_DIRS
  Zoltan_LIBRARY
  Zoltan_LIBRARIES
  )
