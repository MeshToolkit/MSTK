# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE

# Copyright 2019 Triad National Security, LLC. All rights reserved.
#
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE

# -*- mode: cmake -*-
#
# Top-Level CMake file for MSTK
#
# IT IS NOT ADVISABLE TO MODIFY THIS FILE - RATHER ENABLE/DISABLE VARIOUS 
# OPTIONS USING COMMAND LINE ARGUMENTS IN THE config/do-configure-* files
#
# Author: Rao Garimella (rao@lanl.gov)
#

    
project (MSTK C CXX)

cmake_minimum_required(VERSION 3.11)

set(MSTK_VERSION_MAJOR 3)
set(MSTK_VERSION_MINOR 2)
set(MSTK_VERSION_PATCH 0)


cmake_policy(SET CMP0017 NEW)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99")


# Default INSTALL Directory

if (NOT INSTALL_DIR)
  set (INSTALL_DIR ${CMAKE_INSTALL_PREFIX})
endif (NOT INSTALL_DIR)


set (MSTKLIB mstk)


# Add the actual source directory - these will build up the list "mstkfiles"

add_subdirectory(src)


# The main mstk target

add_library(${MSTKLIB} ${mstkfiles})
message(STATUS "Building library " ${MSTKLIB})


# DEBUG FLAG

if (CMAKE_BUILD_TYPE STREQUAL Debug)
  target_compile_definitions(${MSTKLIB} PUBLIC DEBUG)
endif ()


# include files and source files - have this upfront so that an MSTK.h from
# any of the TPL installations doesn't get picked up first

target_include_directories(${MSTKLIB} PUBLIC ${PROJECT_BINARY_DIR}/include)


# If MSTK should enable and use markers to enable fast constant time searches
# in list and set operations

if (MSTK_USE_MARKERS)
  target_compile_definitions(${MSTKLIB} PUBLIC MSTK_USE_MARKERS)
endif ()


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
    set(CMAKE_C_COMPILER ${MPI_C_COMPILER} CACHE FILEPATH "C compiler to use"
      FORCE)
    set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} CACHE FILEPATH "C compiler to use"
      FORCE)
    set(MSTK_HAVE_MPI True CACHE BOOL "If MSTK is compiled with MPI on")
    message(STATUS "MPIEXEC is ${MPIEXEC}")
    message(STATUS "CMAKE_C_COMPILER ${CMAKE_C_COMPILER}")
  else ()
    message(STATUS "Could not find package MPI - build may fail")
    message(STATUS "Try setting CMAKE_C_COMPILER/CMAKE_CXX_COMPILER to the mpi compiler wrappers")
  endif ()
endif ()


################################################################################
# If user wants to read/write ExodusII files, we check if we can locate the
# ExodusII and NetCDF include files and libraries
################################################################################
if (ENABLE_ExodusII)

  target_compile_definitions(${MSTKLIB} PUBLIC ENABLE_ExodusII)

  # For backward compatibility
  if (NOT netCDF_DIR AND NetCDF_DIR)
    set(netCDF_DIR ${NetCDF_DIR})
  endif ()

  ##############################################################################
  # ExodusII - http://sourceforge.net/projects/exodusii/ or 
  # get it as part of Trilinos - http:://trilinos.sandia.gov
  ##############################################################################

  # First seee if a config file got installed as part of the SEACAS project
  # NOTE: NOT ABLE TO PROCESS THIS CORRECTLY
#  find_package(SEACASExodus QUIET CONFIG PATHS ${ExodusII_DIR})

#  if (SEACASExodus_FOUND)
#    target_include_directories(${MSTKLIB} PUBLIC ${ExodusII_INCLUDE_DIRS})
#    target_link_libraries(${MSTKLIB} PUBLIC ${ExodusII_LIBRARIES})
#    message(STATUS "Found Exodus II library (from SEACAS): ${SEACAS_LIBRARIES}")
#  else ()  
    find_package(ExodusII QUIET REQUIRED)  # Uses MSTK module
    
    target_include_directories(${MSTKLIB} PUBLIC ${ExodusII_INCLUDE_DIRS})
    target_link_libraries(${MSTKLIB} PUBLIC ${ExodusII_LIBRARIES})
    message(STATUS "Found Exodus II library: ${ExodusII_LIBRARY}")
#  endif ()

endif ()



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

    find_package(METIS QUIET REQUIRED)   # uses MSTK module
    if (NOT METIS_FOUND)
       message(FATAL_ERROR "Error: could not find Metis library\n"
	  "Try defining Metis_DIR: (-D Metis_DIR:FILEPATH=<Metis_install_prefix>)\n")
    endif ()

    target_include_directories(${MSTKLIB} PUBLIC ${METIS_INCLUDE_DIRS})
    target_compile_definitions(${MSTKLIB} PUBLIC _MSTK_HAVE_METIS)

    if (METIS_MAJOR_VER)
      target_compile_definitions(${MSTKLIB} PUBLIC METIS_${METIS_MAJOR_VER})
    else ()
      target_compile_definitions(${MSTKLIB} PUBLIC METIS_5)
    endif ()

    target_link_libraries(${MSTKLIB} PUBLIC ${METIS_LIBRARIES})

    message(STATUS "Found Metis library: ${METIS_LIBRARY}")
  endif ()

  ##############################################################################
  # Zoltan
  ##############################################################################

  if (ENABLE_ZOLTAN)

    if (ZOLTAN_DIR AND NOT Zoltan_DIR)
      set(Zoltan_DIR ${ZOLTAN_DIR})
    endif()

    # if find_package fails in config mode it will set Zoltan_DIR to
    # Zoltan_DIR-NOTFOUND; so save Zoltan_DIR
    set(zoltan_dir_save ${Zoltan_DIR})
    
    find_package(Zoltan QUIET CONFIG PATHS ${Zoltan_dir})  # try to find thru CMake config file

    if (NOT Zoltan_FOUND)
      set(Zoltan_DIR ${zoltan_dir_save})
      find_package(Zoltan QUIET REQUIRED MODULE)  # fallback to module
    endif ()
  
    target_include_directories(${MSTKLIB} PRIVATE ${Zoltan_INCLUDE_DIRS})
    target_link_libraries(${MSTKLIB} PUBLIC ${Zoltan_LIBRARIES})
    target_compile_definitions(${MSTKLIB} PRIVATE _MSTK_HAVE_ZOLTAN)        
    
    message(STATUS "Found Zoltan library:" ${Zoltan_LIBRARY})
  endif ()
  
endif ()


##############################################################################
# TESTS 
##############################################################################

if (ENABLE_Tests) 

  enable_testing()

  ##############################################################################
  # UnitTest++ - http://unittest-cpp.sourceforge.net/
  ##############################################################################
  
  set(BUILD_TESTS TRUE)

  # For backward compatibility allow UnitTest_DIR for a bit more time
  if (NOT UnitTest++_DIR AND UnitTest_DIR)
    set(UnitTest++_DIR ${UnitTest_DIR})
  endif ()

  set(unittest_dir_save ${UnitTest++_DIR})
  find_package(UnitTest++ QUIET CONFIG PATHS ${UnitTest++_DIR})  # Try to discover thru cmake config file

  if (UnitTest++_FOUND)
    # UnitTest++ sets a weird path for it's includes
    set(UnitTest++_INCLUDE_DIRS ${UTPP_INCLUDE_DIRS}/UnitTest++)

    # Also it sets the target name as UnitTest++ but does not set the
    # LIBRARIES variable
    set(UnitTest++_LIBRARIES UnitTest++)
  else ()
    set(UnitTest++_DIR ${unittest_dir_save})
    find_package(UnitTest++ QUIET REQUIRED MODULE)  # fallback to MSTK module
  endif ()

  message(STATUS "Found UnitTest++: ${UnitTest++_LIBRARY}")

  set(BUILD_TESTS TRUE)
    
endif()


# What are all the TPLs that codes linking to MSTK have to link to

set(MSTK_TPL_LIBRARIES ${METIS_LIBRARIES} ${ZOLTAN_LIBRARIES} 
  ${ExodusII_LIBRARIES} ${NetCDF_LIBRARIES} ${HDF5_LIBRARIES})

# Make a separate variables where the HDF5 libraries are not targets
# but expanded library paths. This will be exported for use by other
# projects

set(MSTK_TPL_LIBRARIES_EXPORT ${METIS_LIBRARIES} ${ZOLTAN_LIBRARIES} 
  ${ExodusII_LIBRARIES} ${NetCDF_LIBRARIES} ${HDF5_LIBRARIES_EXPORT})


#target_link_libraries(${MSTKLIB} PUBLIC ${MSTK_TPL_LIBRARIES})


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


# where to install library
install(TARGETS ${MSTKLIB}
  EXPORT mstk
  ARCHIVE DESTINATION ${INSTALL_DIR}/lib
  LIBRARY DESTINATION ${INSTALL_DIR}/lib
  CONFIGURATIONS ${CMAKE_BUILD_TYPE}
)


#######################################################
# Write a configuration file from template replacing only variables
# enclosed by the @ sign. This will let other programs build on MSTK
# discover how MSTK was built and which TPLs it used. Also put in some
# information into MSTK.h header file
#######################################################

configure_file(${PROJECT_SOURCE_DIR}/cmake/MSTKConfig.cmake.in 
               ${PROJECT_BINARY_DIR}/MSTKConfig.cmake @ONLY)
install(FILES ${PROJECT_BINARY_DIR}/MSTKConfig.cmake 
        DESTINATION ${INSTALL_DIR}/share/cmake/)

file(GLOB ALL_HEADERS ${PROJECT_SOURCE_DIR}/include/*.h)
foreach (HEADER_FILE ${ALL_HEADERS})
  file(COPY ${HEADER_FILE} DESTINATION ${PROJECT_BINARY_DIR}/include)
endforeach(HEADER_FILE)
configure_file(${PROJECT_SOURCE_DIR}/include/MSTK_defines.h 
               ${PROJECT_BINARY_DIR}/include/MSTK_defines.h @ONLY)

if (${INSTALL_DIR} STREQUAL ${PROJECT_SOURCE_DIR})
else ()
INSTALL(FILES 
	${PROJECT_BINARY_DIR}/include/MSTK.h
	${PROJECT_BINARY_DIR}/include/MSTK_defines.h
	${PROJECT_BINARY_DIR}/include/MSTK_types.h
	${PROJECT_BINARY_DIR}/include/MSTK_externs.h
	${PROJECT_BINARY_DIR}/include/MSTK_util.h
        DESTINATION ${INSTALL_DIR}/include)
endif ()



