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
set(MSTK_VERSION_MINOR 83)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")


message("\n-------------------------------------------------------------------")
message("--- Configuring MSTK with:")

################################################################################
# If user wants to read/write ExodusII files, we check if we can locate the
# ExodusII and NetCDF include files and libraries
################################################################################
if (ENABLE_ExodusII)

  add_definitions("-DENABLE_ExodusII")

  ##############################################################################
  # NetCDF - http://www.unidata.ucar.edu/software/netcdf/
  ##############################################################################

# Could replace with find_package(NetCDF REQUIRED)
find_package(NetCDF REQUIRED)
  if(NOT NetCDF_FOUND)
    message(FATAL_ERROR "Error: could not find netcdf library\n"
      "Try defining NetCDF_DIR: (-D NetCDF_DIR:FILEPATH=<NetCDF_install_prefix>)")
  endif()

  message("---- NETCDF_LIBRARY : " ${NetCDF_LIBRARIES})

  ##############################################################################
  # HDF5 - 
  ##############################################################################

 find_package(HDF5 REQUIRED)

 if(NOT HDF5_FOUND)
    message(FATAL_ERROR "Error: could not find hdf5 library\n"
      "Try defining HDF5_DIR: (-D HDF5_DIR:FILEPATH=<HDF5_install_prefix>)")
 endif()


  message("---- HDF5_LIBRARY : " ${HDF5_LIBRARY})

  ##############################################################################
  # ExodusII - http://sourceforge.net/projects/exodusii/ or 
  # get it as part of Trilinos - http:://trilinos.sandia.gov
  ##############################################################################

  find_package(ExodusII REQUIRED)
  if(NOT ExodusII_FOUND)
    message(FATAL_ERROR "Error: could not find ExodusII library\n"
      "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)")
  endif()
  
  include_directories(${ExodusII_INCLUDE_DIRS})
 
  message("---- EXODUSII_LIBRARY : " ${ExodusII_LIBRARIES})

endif(ENABLE_ExodusII)



##############################################################################
# What we need if we need parallel support
##############################################################################

if (ENABLE_PARALLEL)

  if ((NOT ENABLE_METIS) AND (NOT ENABLE_ZOLTAN))
    message(FATAL_ERROR "Error: Must enable METIS or ZOLTAN as a partitioner")
  endif() 

  ##############################################################################
  # Metis
  ##############################################################################

  if (ENABLE_METIS) 

    find_package(METIS REQUIRED)
    if (NOT METIS_FOUND)
       message(FATAL_ERROR "Error: could not find Metis library\n"
	  "Try defining Metis_DIR: (-D Metis_DIR:FILEPATH=<Metis_install_prefix>)\n")

    endif (NOT METIS_FOUND)

    include_directories(${METIS_INCLUDE_DIRS})
  
    message("---- METIS_LIBRARY :" ${METIS_LIBRARY})

    add_definitions(-D_MSTK_HAVE_METIS)

  endif (ENABLE_METIS)

  ##############################################################################
  # Zoltan
  ##############################################################################

  if (ENABLE_ZOLTAN) 

    find_package(ZOLTAN REQUIRED)
    if (NOT ZOLTAN_FOUND)
       message(FATAL_ERROR "Error: could not find Zoltan library\n"
	  "Try defining Zoltan_DIR: (-D Zoltan_DIR:FILEPATH=<Zoltan_install_prefix>)\n")
     endif (NOT ZOLTAN_FOUND)
  
     include_directories(${ZOLTAN_INCLUDE_DIRS})
  
     message("---- ZOLTAN_LIBRARY :" ${ZOLTAN_LIBRARY})

     add_definitions(-D_MSTK_HAVE_ZOLTAN)

  endif (ENABLE_ZOLTAN)

endif (ENABLE_PARALLEL) 


##############################################################################
# TESTS 
##############################################################################

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
      HINTS ${UnitTest_DIR}
      ${UnitTest_LIB_DIR}
      ${UnitTest_DIR}/lib
      )
    if (NOT UNITTEST_LIBRARY)
      message(FATAL_ERROR "Error: could not find unittest++ library")
    endif()
  
    find_path(UNITTEST_INCLUDES "UnitTest++.h" 
      HINTS ${UnitTest_INCLUDE_DIR}
      ${UnitTest_DIR}
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


#############################################################################
# Now add the subdirectory containing the source
#############################################################################

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
message(STATUS "Building library " ${MSTKLIB})


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
# Unittests
###############################################################################

if (ENABLE_Tests)
  add_subdirectory(unittests)
endif ()

#############################################################################
# Utilities
###############################################################################

add_subdirectory(utils)


#############################################################################
# Example
###############################################################################

add_subdirectory(example)


