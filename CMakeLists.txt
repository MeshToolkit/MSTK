# -*- mode: cmake -*-
#
# Top-Level CMake file for MSTK
#
# IT IS NOT ADVISABLE TO MODIFY THIS FILE - RATHER ENABLE/DISABLE VARIOUS 
# OPTIONS USING COMMAND LINE ARGUMENTS IN THE config/do-configure-* files
#
# Author: Rao Garimella (rao@lanl.gov)
#

    
project (MSTK)

cmake_minimum_required(VERSION 2.8)

set(MSTK_VERSION_MAJOR 1)
set(MSTK_VERSION_MINOR 8)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})


################################################################################
# If user wants to read/write ExodusII files, we check if we can locate the
# ExodusII and NetCDF include files and libraries
################################################################################

if (ENABLE_ExodusII)

  add_definition("-DENABLE_ExodusII")

  ##############################################################################
  # NetCDF - http://www.unidata.ucar.edu/software/netcdf/
  ##############################################################################

  find_library(NETCDF_LIBRARY netcdf
    HINTS ${NetCDF_LIB_DIR}
    ${NetCDF_DIR}/lib
    /usr/lib
    /usr/local/lib
    )

  if(NOT NETCDF_LIBRARY)
    message(FATAL_ERROR "Error: could not find netcdf library\n"
      "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
      "\nor defining NetCDF_LIB_DIR")
  endif()

  find_path(NETCDF_INCLUDES netcdf.h 
    HINTS ${NetCDF_INCLUDE_DIR}
    ${NetCDF_DIR}/include
    /usr/include
    /usr/local/include
    )
  if(NOT NETCDF_INCLUDES)
    message(FATAL_ERROR "Error: could not find NetCDF includes\n"
      "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)"
      "\nor defining NetCDF_INCLUDE_DIR")
  endif()


  include_directories(${NETCDF_INCLUDES})

  #
  # DID NOT SEEM TO NEED THIS - but left uncommented in checked in version
  #
  add_library(netcdf SHARED IMPORTED)
  set_property(TARGET netcdf PROPERTY IMPORTED_LOCATION ${NETCDF_LIBRARY})

  #target_link_libraries(netcdf )   # See note under ExodusII.
  set(NETCDF_LIBS_ALL netcdf hdf5 hdf5_hl z)


  ##############################################################################
  # ExodusII - http://sourceforge.net/projects/exodusii/ or 
  # get it as part of Trilinos - http:://trilinos.sandia.gov
  ##############################################################################

  find_library(EXODUSII_LIBRARY exoIIv2c
    HINTS ${ExodusII_LIB_DIR}
    ${ExodusII_DIR}/cbind/src
    ${ExodusII_DIR}/lib
    )

  if(NOT EXODUSII_LIBRARY)
    message(FATAL_ERROR "Error: could not find ExodusII library\n"
      "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)\n"
      "Or ExodusII_LIB_DIR")
  endif()
  
  find_path(EXODUSII_INCLUDES exodusII.h 
    HINTS ${ExodusII_INCLUDE_DIR}
    ${ExodusII_DIR}/cbind/include
    ${ExodusII_DIR}/include
    )
  
  if(NOT EXODUSII_INCLUDES)
    message(FATAL_ERROR "Error: could not find ExodusII includes\n"
      "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)\n"
      "Or ExodusII_INCLUDE_DIR")
  endif()

  include_directories(${EXODUSII_INCLUDES})
 
  add_library(exodusii STATIC IMPORTED)
  set_property(TARGET exodusii PROPERTY IMPORTED_LOCATION ${EXODUSII_LIBRARY})
  
  # Can't add link dependencies to imported libraries. See http://www.vtk.org/Bug/view.php?id=10395
  # target_link_libraries(exodusii netcdf)
  
  # Instead, we create a list of all the libraries required for ExodusII to link.
  set(EXODUS_LIBS_ALL exodusii netcdf) 

endif(ENABLE_ExodusII)



##############################################################################
# What we need if we need parallel support
##############################################################################

if (ENABLE_PARALLEL)

  ##############################################################################
  # Metis
  ##############################################################################

  find_library(METIS_LIBRARY metis
    HINTS ${METIS_LIB_DIR}
    ${Metis_DIR}
    ${Metis_DIR}/Lib
    )

  if(NOT METIS_LIBRARY)
    message(FATAL_ERROR "Error: could not find Metis library\n"
      "Try defining Metis_DIR: (-D Metis_DIR:FILEPATH=<Metis_install_prefix>)\n"
      "Or METIS_LIB_DIR")
  endif()
  
  find_path(METIS_INCLUDES metis.h 
    HINTS ${METIS_INCLUDE_DIR}
    ${Metis_DIR}
    ${Metis_DIR}/include
    ${Metis_DIR}/Lib
    )
  
  if(NOT METIS_INCLUDES)
    message(FATAL_ERROR "Error: could not find Metis includes\n"
      "Try defining Metis_DIR: (-D Metis_DIR:FILEPATH=<Metis_install_prefix>)\n"
      "Or Metis_INCLUDE_DIR")
  endif()

  include_directories(${METIS_INCLUDES})
  
  add_library(metis STATIC IMPORTED)
  set_property(TARGET metis PROPERTY IMPORTED_LOCATION ${METIS_LIBRARY})
  
  set(METIS_LIBS_ALL metis) 

endif (ENABLE_PARALLEL)



if (ENABLE_Tests) 

  enable_testing()

  ##############################################################################
  # UnitTest++ - http://unittest-cpp.sourceforge.net/
  ##############################################################################
  
  if (NOT DEFINED UnitTest_DIR)
    message(WARNING "No UnitTest++ directory defined. To compile the unit tests add:
 -D UnitTest_DIR::FILEPATH=<UnitTest++_install_prefix>")
  else()
    find_library(UNITTEST_LIBRARY 
      NAMES UnitTest++ unittest++
      HINTS ${UnitTest_LIB_DIR}
      ${UnitTest_DIR}/lib
      )
    if (NOT UNITTEST_LIBRARY)
      message(FATAL_ERROR "Error: could not find unittest++ library")
    endif()
  
    find_path(UNITTEST_INCLUDES "UnitTest++.h" 
      HINTS ${UnitTest_INCLUDE_DIR}
      ${UnitTest_DIR}/include
      ${UnitTest_DIR}/src 
      ${UnitTest_DIR}/include/unittest++
      )
    if (NOT UNITTEST_INCLUDES)
      message(FATAL_ERROR "Error: could not find UnitTest++ includes")
    else()
      include_directories(${UNITTEST_INCLUDES})
    endif()
    
    set(BUILD_TESTS 1)
    add_library(unittest++ STATIC IMPORTED)
    set_property(TARGET unittest++ PROPERTY IMPORTED_LOCATION ${UNITTEST_LIBRARY})

  endif()
endif()




##############################################################################
# MAIN MSTK LIBRARY BUILD
##############################################################################

# Flags for parallel build

if (ENABLE_PARALLEL)
  add_definitions(-DMSTK_HAVE_MPI)
endif (ENABLE_PARALLEL)


# Default build type is Release (Optimized)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()


# Default INSTALL Directory

if (NOT INSTALL_DIR)
  set (INSTALL_DIR ${PROJECT_SOURCE_DIR})
endif ()


# include files and source files

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)


# Extra library suffix for Debug builds

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  set (CONFIG_SUFFIX -d)
endif ()

# Add version number to library name

if (INSTALL_ADD_VERSION)
  set (MSTKLIB mstk-${MSTK_VERSION_MAJOR}.${MSTK_VERSION_MINOR}${CONFIG_SUFFIX})
else ()
  set (MSTKLIB mstk${CONFIG_SUFFIX})
endif ()

add_library(${MSTKLIB} ${mstkfiles})


# where to install library

install(TARGETS ${MSTKLIB}
  ARCHIVE
  DESTINATION ${INSTALL_DIR}/lib/${ARCHOS}
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
  )

# export cmake file other projects can use to find MSTK

export(TARGETS ${MSTKLIB} FILE ${INSTALL_DIR}/mstk.cmake)
message(STATUS ${MSTKLIB})


# where to install include files

if (${INSTALL_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
else ()
INSTALL(FILES 
	${PROJECT_SOURCE_DIR}/include/MSTK.h
	${PROJECT_SOURCE_DIR}/include/MSTK_defines.h
	${PROJECT_SOURCE_DIR}/include/MSTK_types.h
	${PROJECT_SOURCE_DIR}/include/MSTK_externs.h
	${PROJECT_SOURCE_DIR}/include/MSTK_util.h
	${PROJECT_SOURCE_DIR}/include/MSTK_malloc.h
        DESTINATION ${INSTALL_DIR}/include)
endif ()

# Don't know how to export the include files locations to mstk.cmake


#############################################################################
# Example
###############################################################################

add_subdirectory(example)
