# Copyright: 2019- Triad National Security, LLC
#
# ExodusII Find Module for MSTK
#
# ExodusII needs METIS; This module will try to find METIS as well and add it
# as a dependency
#
# Usage: To search a particular path you can specify the path in
# CMAKE_PREFIX_PATH, in the CMake variable ExodusII_DIR or environment
# variable ExodusII_ROOT
#
# Following variables are set:
# ExodusII_FOUND          (BOOL)   Flag indicating if ExodusII was found
# ExodusII_INCLUDE_DIRS   (PATH)   Path to ExodusII include files
# ExodusII_LIBRARY        (FILE)   ExodusII library (libzoltan.a, libzoltan.so)
# ExodusII_LIBRARIES      (LIST)   List of ExodusII targets (ExodusII::ExodusII)
# ExodusII_ROOT           (PATH)   Top level directory where Exodus is installed
#
#
# Additional variables
# ExodusII_VERSION          (STRING)     ExodusII Version string
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_ExodusII QUIET exodusII)


# Search for include files
find_path(ExodusII_INCLUDE_DIR
  NAMES exodusII.h
  HINTS ${PC_ExodusII_INCLUDE_DIRS}
  PATHS ${ExodusII_DIR}
  PATH_SUFFIXES include
  )

if (NOT ExodusII_INCLUDE_DIR)
  if (ExodusII_FIND_REQUIRED)
    message(FATAL "Cannot locate exodusII.h")
  else ()
    if (NOT ExodusII_FIND_QUIETLY)
      message(WARNING "Cannot locate exodusII.h")
    endif ()
  endif ()
endif ()

set(ExodusII_INCLUDE_DIRS "${ExodusII_INCLUDE_DIR}")

# Search for libraries

find_library(ExodusII_LIBRARY
  NAMES exodus exoIIv2c
  HINTS ${PC_ExodusII_LIBRARY_DIRS}
  PATHS ${ExodusII_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT ExodusII_LIBRARY)
  if (ExodusII_FIND_REQUIRED)
    message(FATAL "Can not locate ExodusII library")
  else (ExodusII_FIND_REQUIRED)
    if (NOT ExodusII_FIND_QUIETLY)
      message(WARNING "Cannot locate ExodusII library")
    endif ()
  endif ()
endif ()

# Set library version

set(ExodusII_VERSION ${PC_ExodusII_VERSION})  # No guarantee
if (NOT ExodusII_VERSION AND ExodusII_INCLUDE_DIR)
  set(exodus_h "${ExodusII_INCLUDE_DIR}/exodusII.h")
  file(STRINGS "${exodus_h}" exodus_version_string REGEX "^#define EX_API_VERS")
  string(REGEX REPLACE "^#define EX_API_VERS ([0-9]+\\.[0-9]+).*$" "\\1" exodus_version "${exodus_version_string}")

  set(ExodusII_VERSION "${exodus_version}")
endif ()

# Not sure if this is the right way to do it, but this is to help
# other upstream packages that attempt to find the ExodusII package
# due to transitive dependencies
if (NOT ExodusII_ROOT)
  get_filename_component(ExodusII_ROOT "${ExodusII_INCLUDE_DIR}/.." ABSOLUTE CACHE "Top level dir of ExodusII installation" FORCE)
endif ()
if (NOT ExodusII_DIR)
  get_filename_component(ExodusII_DIR "${ExodusII_INCLUDE_DIR}/.." ABSOLUTE CACHE "Top level dir of ExodusII installation" FORCE)
endif ()


# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(ExodusII
  DEFAULT_MSG
  ExodusII_LIBRARY ExodusII_INCLUDE_DIR ExodusII_ROOT)

# find_package_handle_standard_args ignores case and sets PACKAGE_FOUND
if (NOT ExodusII_FOUND AND EXODUSII_FOUND)
  set(ExodusII_FOUND ${EXODUSII_FOUND})
endif ()

# Create ExodusII target

if (ExodusII_FOUND AND NOT TARGET ExodusII::ExodusII)
  set(ExodusII_LIBRARIES ExodusII::ExodusII)
  add_library(${ExodusII_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${ExodusII_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${ExodusII_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_ExodusII_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${ExodusII_INCLUDE_DIR}")

  # ExodusII depends on netCDF. Attempt to find it
  
  # Try to discover through cmake config file (usually named
  # netcdf-config.cmake netCDFConfig.cmake)
  find_package(netCDF NAMES netcdf netCDF CONFIG HINTS ${netCDF_DIR})
  if (NOT netCDF_FOUND)
    # Fallback to MSTK module named FindNetCDF.cmake
    find_package(netCDF QUIET REQUIRED MODULE)
  endif ()

  # If the package was found using the config file, only netCDF_DIR
  # may be set not netCDF_ROOT.
  if (netCDF_DIR AND NOT netCDF_ROOT)
    set(netCDF_ROOT ${netCDF_DIR} CACHE PATH "Top level installation dir of netCDF")
  endif ()

  # Add netCDF as a dependency of ExodusII
  target_link_libraries(${ExodusII_LIBRARIES} INTERFACE ${netCDF_LIBRARIES})
endif()

  
# Hide these variables from the cache
mark_as_advanced(
  ExodusII_INCLUDE_DIR
  ExodusII_INCLUDE_DIRS
  ExodusII_LIBRARY
  ExodusII_LIBRARIES
  ExodusII_VERSION
  )

