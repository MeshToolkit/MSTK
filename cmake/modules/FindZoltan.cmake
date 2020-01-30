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
# Zoltan_LIBRARIES      (LIST)   List of Zoltan targets (Zoltan::Zoltan)
# Zoltan_ROOT           (PATH)   Top level directory where Zoltan is installed
# Zoltan_DIR            (PATH)   Top level directory where Zoltan is installed
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_Zoltan QUIET zoltan)


# Search for include files

find_path(Zoltan_INCLUDE_DIR
  NAMES zoltan.h
  HINTS ${PC_Zoltan_INCLUDE_DIRS}
  PATHS ${Zoltan_DIR}
  PATH_SUFFIXES include)

if (NOT Zoltan_INCLUDE_DIR)
  if (Zoltan_FIND_REQUIRED)
    message(FATAL "Cannot locate zoltan.h")
  else (Zoltan_FIND_REQUIRED)
    if (NOT Zoltan_FIND_QUIETLY)
      message(WARNING "Cannot locate zoltan.h")
    endif ()
  endif ()
endif ()

set(Zoltan_INCLUDE_DIRS "${Zoltan_INCLUDE_DIR}")


# Search for libraries

find_library(Zoltan_LIBRARY
  NAMES zoltan
  HINTS ${PC_Zoltan_LIBRARY_DIRS}
  PATHS ${Zoltan_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT Zoltan_LIBRARY)
  if (Zoltan_FIND_REQUIRED)
    message(FATAL "Can not locate Zoltan library")
  else (Zoltan_FIND_REQUIRED)
    if (NOT Zoltan_FIND_QUIETLY)
      message(WARNING "Cannot locate Zoltan library")
    endif ()
  endif ()
endif ()

set(Zoltan_VERSION PC_Zoltan_VERSION})  # No guarantee

# Not sure if this is the right way to do it, but this is to help
# other upstream packages that attempt to find the Zoltan package
# due to transitive dependencies
if (NOT Zoltan_ROOT)
  set(Zoltan_ROOT "${Zoltan_INCLUDE_DIR}/.." CACHE PATH "Top level dir of Zoltan installation" FORCE)
endif ()


# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Zoltan
  DEFAULT_MSG
  Zoltan_LIBRARY Zoltan_INCLUDE_DIR Zoltan_ROOT)

# find_package_handle_standard_args ignores case and sets PACKAGE_FOUND
if (NOT Zoltan_FOUND AND ZOLTAN_FOUND)
  set(Zoltan_FOUND ${ZOLTAN_FOUND})
endif ()

# Create Zoltan target

if (Zoltan_FOUND AND NOT TARGET Zoltan::Zoltan)
  set(Zoltan_LIBRARIES Zoltan::Zoltan)
  add_library(${Zoltan_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${Zoltan_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${Zoltan_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_Zoltan_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${Zoltan_INCLUDE_DIR}")

  # Don't know how to automatically determine if Zoltan needs ParMETIS
  # or PTScotch. So, try to find these packages and add them as
  # dependencies if they are found. Cannot do any harm
  
  
  # ParMETIS does not (and likely will not) have a CMake config file
  # but one can dream
  
  set(parmetis_dir_save ${ParMETIS_DIR})
  find_package(ParMETIS QUIET CONFIG)
  if (NOT ParMETIS_FOUND)
    set(ParMETIS_DIR ${parmetis_dir_save})
    find_package(ParMETIS QUIET MODULE)
  endif ()
  
  if (ParMETIS_FOUND)
    # Add ParMETIS as a dependency of Zoltan
    target_link_libraries(${Zoltan_LIBRARIES} INTERFACE ${ParMETIS_LIBRARIES})
  endif ()
  
  
  # PTScotch does not have a CMake config file (one cannot even dream)
  
  set(ptscotch_dir_save ${PTScotch_DIR})
  find_package(PTScotch QUIET CONFIG)
  if (NOT PTScotch_FOUND)
    set(PTScotch_DIR ${ptscotch_dir_save})
    find_package(PTScotch QUIET MODULE)
  endif ()
  
  if (PTScotch_FOUND)
    # Add PTScotch as a dependency of Zoltan
    target_link_libraries(${Zoltan_LIBRARIES} INTERFACE ${PTScotch_LIBRARIES})
  endif ()
endif ()

# Hide these variables from the cache
mark_as_advanced(
  Zoltan_INCLUDE_DIR
  Zoltan_INCLUDE_DIRS
  Zoltan_LIBRARY
  Zoltan_LIBRARIES
  )
