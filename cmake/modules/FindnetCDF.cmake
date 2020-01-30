# Copyright: 2019- Triad National Security, LLC
#
# netCDF Find Module for MSTK
#
# Usage: To search a particular path you can specify the path in
# CMAKE_PREFIX_PATH, in the CMake variable netCDF_DIR or environment
# variable netCDF_ROOT
#
# Following variables are set:
# netCDF_FOUND          (BOOL)   Flag indicating if netCDF was found
# netCDF_INCLUDE_DIRS   (PATH)   Path to netCDF include files
# netCDF_LIBRARY        (FILE)   netCDF library (libnetcdf.a, libnetcdf.so)
# netCDF_LIBRARIES      (LIST)   List of netCDF targets (netCDF::netCDF)
# netCDF_ROOT           (PATH)   Top level directory where netCDF is installed
#
# #############################################################################

# First use pkg-config to parse an installed .pc file to find the
# library although we cannot rely on it

find_package(PkgConfig)
pkg_check_modules(PC_netCDF QUIET netCDF)


# Search for include files

find_path(netCDF_INCLUDE_DIR
  NAMES netcdf.h
  HINTS ${PC_netCDF_INCLUDE_DIRS} ${netCDF_DIR}
  PATH_SUFFIXES include)

if (NOT netCDF_INCLUDE_DIR)
  if (netCDF_FIND_REQUIRED)
    message(FATAL "Cannot locate netcdf.h")
  else (netCDF_FIND_REQUIRED)
    if (NOT netCDF_FIND_QUIETLY)
      message(WARNING "Cannot locate netcdf.h")
    endif ()
  endif ()
endif ()

set(netCDF_INCLUDE_DIRS "${netCDF_INCLUDE_DIR}")


# Search for libraries

find_library(netCDF_LIBRARY
  NAMES netcdf
  HINTS ${PC_netCDF_LIBRARY_DIRS} ${netCDF_DIR}
  PATH_SUFFIXES lib lib64)

if (NOT netCDF_LIBRARY)
  if (netCDF_FIND_REQUIRED)
    message(FATAL "Can not locate netCDF library")
  else (netCDF_FIND_REQUIRED)
    if (NOT netCDF_FIND_QUIETLY)
      message(WARNING "Cannot locate netCDF library")
    endif ()
  endif ()
endif ()

set(netCDF_VERSION PC_netCDF_VERSION})  # No guarantee

# Not sure if this is the right way to do it, but this is to help
# other upstream packages that attempt to find the netCDF package
# due to transitive dependencies
if (NOT netCDF_ROOT)
  get_filename_component(netCDF_ROOT "${netCDF_INCLUDE_DIR}/.." ABSOLUTE CACHE "Top level dir of netCDF installation" FORCE)
endif ()
if (NOT netCDF_DIR)
  get_filename_component(netCDF_DIR "${netCDF_INCLUDE_DIR}/.." ABSOLUTE CACHE "Top level dir of netCDF installation" FORCE)
endif ()


# Finish setting standard variables if everything is found
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(netCDF
  FOUND_VAR netCDF_FOUND
  REQUIRED_VARS netCDF_LIBRARY netCDF_INCLUDE_DIR netCDF_ROOT)

# find_package_handle_standard_args ignores case and sets PACKAGE_FOUND
if (NOT netCDF_FOUND AND NETCDF_FOUND)
  set(netCDF_FOUND ${NETCDF_FOUND})
endif ()

# Create netCDF target

if (netCDF_FOUND AND NOT TARGET netCDF::netCDF)
  set(netCDF_LIBRARIES netCDF::netCDF)
  add_library(${netCDF_LIBRARIES} UNKNOWN IMPORTED)
  set_target_properties(${netCDF_LIBRARIES} PROPERTIES
    IMPORTED_LOCATION "${netCDF_LIBRARY}"
    INTERFACE_COMPILE_OPTIONS "${PC_netCDF_CFLAGS_OTHER}"
    INTERFACE_INCLUDE_DIRECTORIES "${netCDF_INCLUDE_DIR}")

  # netCDF depends on HDF5
  
  find_package(HDF5 QUIET NAMES hdf5 COMPONENTS C HL CONFIG)
  
  if (NOT HDF5_FOUND)
    find_package(HDF5 QUIET REQUIRED COMPONENTS C HL MODULE)  # Fallback to built in module
  endif ()

  target_link_libraries(${netCDF_LIBRARIES} INTERFACE ${HDF5_LIBRARIES})
endif ()


# Hide these variables from the cache
mark_as_advanced(
  netCDF_INCLUDE_DIR
  netCDF_INCLUDE_DIRS
  netCDF_LIBRARY
  netCDF_LIBRARIES
  )
