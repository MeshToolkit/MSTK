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

cmake_policy(SET CMP0017 NEW)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")

# This option may be needed for machines like Edison and Hopper

if (PREFER_STATIC_LIBRARIES)
   set(CMAKE_FIND_LIBRARY_SUFFIXES ".a; .so")
endif()

# Default build type is Release (Optimized)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release)
endif ()


###############################################################################
# include files and source files - have this upfront so that an MSTK.h from
# any of the TPL installations doesn't get picked up first
##############################################################################

include_directories(${PROJECT_SOURCE_DIR}/include)

add_subdirectory(src)
################################################################################
# Discover TPLs (Third Party Libraries) based on configure options
################################################################################

message("\n-------------------------------------------------------------------")
message("--- Configuring MSTK with:")

#
# Flags for parallel build

if (ENABLE_PARALLEL)
  find_package(MPI)
  if (MPI_FOUND)
    message("PARALLELISM ON")
    add_definitions(-DMSTK_HAVE_MPI)
    message(STATUS "MPIEXEC is ${MPIEXEC}")
  else (MPI_FOUND)
    message("PARALLELISM OFF")
  endif (MPI_FOUND)
endif (ENABLE_PARALLEL)


################################################################################
# If MSTK should enable and use markers to enable fast constant time searches
# in list and set operations
################################################################################

if (MSTK_USE_MARKERS)
  add_definitions("-DMSTK_USE_MARKERS")
endif (MSTK_USE_MARKERS)

################################################################################
# If user wants to read/write ExodusII files, we check if we can locate the
# ExodusII and NetCDF include files and libraries
################################################################################
if (ENABLE_ExodusII)

  add_definitions("-DENABLE_ExodusII")

  if (EXODUS_6_DEPRECATED)
    add_definitions("-DEXODUS_6_DEPRECATED")
  endif (EXODUS_6_DEPRECATED)

  ##############################################################################
  # HDF5 - 
  ##############################################################################

 set(HDF5_ROOT ${HDF5_DIR})
 find_package(HDF5 REQUIRED 
              COMPONENTS C CXX HL)

 if(NOT HDF5_FOUND)
    message(FATAL_ERROR "Error: could not find hdf5 library\n"
      "Try defining HDF5_DIR: (-D HDF5_DIR:FILEPATH=<HDF5_install_prefix>)")
 endif()


  message("---- HDF5_LIBRARIES : " ${HDF5_LIBRARIES})
  message("---- HDF5_IS_PARALLEL : " ${HDF5_IS_PARALLEL})

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

    if (METIS_MAJOR_VER)
      add_definitions(-DMETIS_${METIS_MAJOR_VER})
    else (METIS_MAJOR_VER)
      add_definitions(-DMETIS_5)
    endif (METIS_MAJOR_VER)

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

     if (ZOLTAN_NEEDS_ParMETIS)
       find_package(ParMETIS REQUIRED)
       if (NOT ParMETIS_FOUND)
         message(FATAL_ERROR  "Error: could not find ParMETIS library \
	    needed by ZOLTAN\n \
	   Try defining ParMETIS_DIR: \
	   (-D ParMETIS_DIR:FILEPATH=<ParMETIS_install_prefix>)\n")
       endif (NOT ParMETIS_FOUND)

       list(APPEND ZOLTAN_LIBRARIES ${ParMETIS_LIBRARIES})
       message("ParMETIS_LIBRARIES --- " ${ParMETIS_LIBRARIES})
       message("ZOLTAN_LIBRARIES --- " ${ZOLTAN_LIBRARIES})
     endif (ZOLTAN_NEEDS_ParMETIS)

     if (ZOLTAN_NEEDS_PTScotch)
       find_package(PTScotch REQUIRED)
       if (NOT PTScotch_FOUND)
         message(FATAL_ERROR  "Error: could not find PTScotch library \
	    needed by ZOLTAN\n \
	   Try defining PTScotch_DIR: \
	   (-D PTScotch_DIR:FILEPATH=<PTScotch_install_prefix>)\n")
       endif (NOT PTScotch_FOUND)

       list(APPEND ZOLTAN_LIBRARIES ${PTScotch_LIBRARIES})
       message("PTScotch_LIBRARIES --- " ${PTScotch_LIBRARIES})
       message("ZOLTAN_LIBRARIES --- " ${ZOLTAN_LIBRARIES})
     endif (ZOLTAN_NEEDS_PTScotch)

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
      ${UnitTest_DIR}/include/UnitTest++
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


# What are all the TPLs that codes linking to MSTK have to link to

set(MSTK_TPL_LIBRARIES ${METIS_LIBRARIES} ${ZOLTAN_LIBRARIES} 
  ${ExodusII_LIBRARIES} ${NetCDF_LIBRARIES} ${HDF5_LIBRARIES})

# Make a separate variables where the HDF5 libraries are not targets
# but expanded library paths. This will be exported for use by other
# projects

set(MSTK_TPL_LIBRARIES_EXPORT ${METIS_LIBRARIES} ${ZOLTAN_LIBRARIES} 
  ${ExodusII_LIBRARIES} ${NetCDF_LIBRARIES} ${HDF5_LIBRARIES_EXPORT})


##############################################################################
# MSTK LIBRARY BUILD CONFIGURATION
##############################################################################


# Extra library suffix for Debug builds

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  if (INSTALL_ADD_DEBUG_SUFFIX)
    set (CONFIG_SUFFIX -d)
  endif ()
  add_definitions(-DDEBUG)
endif ()

# Add version number to library name

if (INSTALL_ADD_VERSION)
  set (MSTKLIB mstk-${MSTK_VERSION_MAJOR}.${MSTK_VERSION_MINOR}${CONFIG_SUFFIX})
else ()
  set (MSTKLIB mstk${CONFIG_SUFFIX})
endif ()

add_library(${MSTKLIB} ${mstkfiles})
message(STATUS "Building library " ${MSTKLIB})

# Useful to define this particularly for shared libraries on OSX

target_link_libraries(${MSTKLIB} ${MSTK_TPL_LIBRARIES})


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


##########################################################
# Configuration of installation options
##########################################################

# Default INSTALL Directory

if (NOT INSTALL_DIR)
  set (INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
endif ()

# where to install library
if (INSTALL_PREFIX_ARCHOS) 
  install(TARGETS ${MSTKLIB}
    EXPORT mstk
    ARCHIVE DESTINATION ${INSTALL_DIR}/lib/${ARCHOS}
    LIBRARY DESTINATION ${INSTALL_DIR}/lib/${ARCHOS}
    CONFIGURATIONS ${CMAKE_BUILD_TYPE}
  )
else ()
  install(TARGETS ${MSTKLIB}
    EXPORT mstk
    ARCHIVE DESTINATION ${INSTALL_DIR}/lib
    LIBRARY DESTINATION ${INSTALL_DIR}/lib
    CONFIGURATIONS ${CMAKE_BUILD_TYPE}
  )
endif ()

# where to install include files

if (${INSTALL_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
else ()
INSTALL(FILES 
	${PROJECT_SOURCE_DIR}/include/MSTK.h
	${PROJECT_SOURCE_DIR}/include/MSTK_defines.h
	${PROJECT_SOURCE_DIR}/include/MSTK_types.h
	${PROJECT_SOURCE_DIR}/include/MSTK_externs.h
	${PROJECT_SOURCE_DIR}/include/MSTK_util.h
        DESTINATION ${INSTALL_DIR}/include)
endif ()

############################################################################## 
# Write a configuration file from template replacing only variables enclosed
# by the @ sign. This will let other programs build on MSTK discover how
# MSTK was built and which TPLs it used
#############################################################################

configure_file(${PROJECT_SOURCE_DIR}/cmake/MSTKConfig.cmake.in 
               ${PROJECT_BINARY_DIR}/MSTKConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/MSTKConfig.cmake 
        DESTINATION ${INSTALL_DIR}/share/cmake/)


