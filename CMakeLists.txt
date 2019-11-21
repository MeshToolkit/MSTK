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
endif (MSTK_USE_MARKERS)


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
  else (MPI_FOUND)
    message(STATUS "Could not find package MPI - build may fail")
    message(STATUS "Try setting CMAKE_C_COMPILER/CMAKE_CXX_COMPILER to the mpi compiler wrappers")
  endif (MPI_FOUND)
endif (ENABLE_PARALLEL)


################################################################################
# If user wants to read/write ExodusII files, we check if we can locate the
# ExodusII and NetCDF include files and libraries
################################################################################
if (ENABLE_ExodusII)

  target_compile_definitions(${MSTKLIB} PUBLIC ENABLE_ExodusII)

  # For backward compatibility
  if (NOT netCDF_DIR AND NetCDF_DIR)
    set(netCDF_DIR ${NetCDF_DIR})
  endif (NOT netCDF_DIR AND NetCDF_DIR)

  ##############################################################################
  # ExodusII - http://sourceforge.net/projects/exodusii/ or 
  # get it as part of Trilinos - http:://trilinos.sandia.gov
  ##############################################################################

  find_package(ExodusII REQUIRED)  # Uses MSTK module
  if(NOT ExodusII_FOUND)
    message(FATAL_ERROR "Error: could not find ExodusII library\n"
      "Try defining ExodusII_DIR: (-D ExodusII_DIR:FILEPATH=<ExodusII_install_prefix>)")
  endif(NOT ExodusII_FOUND)
  
  target_include_directories(${MSTKLIB} PUBLIC ${ExodusII_INCLUDE_DIRS})
  target_link_libraries(${MSTKLIB} PUBLIC ${ExodusII_LIBRARIES})

  message(STATUS "### Exodus II library: ${ExodusII_LIBRARIES}")

  get_target_property(ExodusII_LOC ${ExodusII_LIBRARIES} IMPORTED_LOCATION)
  if (ExodusII_LOC)
    message(STATUS "### Exodus II location: ${ExodusII_LOC}")
  else (ExodusII_LOC)
    if (ExodusII_LIBRARY_DIRS)
      message(STATUS "### Exodus II library dir: ${ExodusII_LIBRARY_DIRS}")
    endif (ExodusII_LIBRARY_DIRS)
  endif (ExodusII_LOC)
 
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

    find_package(METIS REQUIRED)   # uses MSTK module
    if (NOT METIS_FOUND)
       message(FATAL_ERROR "Error: could not find Metis library\n"
	  "Try defining Metis_DIR: (-D Metis_DIR:FILEPATH=<Metis_install_prefix>)\n")

    endif (NOT METIS_FOUND)

    target_include_directories(${MSTKLIB} PUBLIC ${METIS_INCLUDE_DIRS})
    target_compile_definitions(${MSTKLIB} PUBLIC _MSTK_HAVE_METIS)

    if (METIS_MAJOR_VER)
      target_compile_definitions(${MSTKLIB} PUBLIC METIS_${METIS_MAJOR_VER})
    else (METIS_MAJOR_VER)
      target_compile_definitions(${MSTKLIB} PUBLIC METIS_5)
    endif (METIS_MAJOR_VER)

    target_link_libraries(${MSTKLIB} PUBLIC ${METIS_LIBRARIES})

    message(STATUS "### Metis library: ${METIS_LIBRARIES}")
    get_target_property(METIS_LOC ${METIS_LIBRARIES} IMPORTED_LOCATION)
    if (METIS_LOC)
      message(STATUS "### METIS location: ${METIS_LOC}")
    else (METIS_LOC)
      if (METIS_LIBRARY_DIRS)
	message(STATUS "### METIS library dir: ${METIS_LIBRARY_DIRS}")
      endif (METIS_LIBRARY_DIRS)
    endif (METIS_LOC)
  endif (ENABLE_METIS)

  ##############################################################################
  # Zoltan
  ##############################################################################

  if (ENABLE_ZOLTAN) 

    find_package(Zoltan REQUIRED HINTS ${ZOLTAN_DIR})
    if (Zoltan-NOTFOUND)
      message(FATAL_ERROR "Error: could not find Zoltan library\n"
	"Try defining Zoltan_DIR: (-D Zoltan_DIR:FILEPATH=<Zoltan_install_prefix>)\n")
    endif (Zoltan-NOTFOUND)
  

    target_include_directories(${MSTKLIB} PRIVATE ${Zoltan_INCLUDE_DIRS})
    
    message(STATUS "### Zoltan library:" ${Zoltan_LIBRARIES})
    get_target_property(Zoltan_LOC ${Zoltan_LIBRARIES} IMPORTED_LOCATION)
    if (Zoltan_LOC)
      message(STATUS "### Zoltan location: ${Zoltan_LOC}")
    else (Zoltan_LOC)
      if (Zoltan_LIBRARY_DIRS)
	message(STATUS "### Zoltan library dir: ${Zoltan_LIBRARY_DIRS}")
      endif (Zoltan_LIBRARY_DIRS)
    endif (Zoltan_LOC)
    
    target_compile_definitions(${MSTKLIB} PUBLIC MSTK_HAVE_ZOLTAN)

    target_link_libraries(${MSTKLIB} PUBLIC ${Zoltan_LIBRARIES})
    
    # Code we introduced for MSTK to work in XSDK
    # If we use targets we shouldn't have to do this manually
    # BUT KEEP COMMENTED OUT JUST IN CASE WE HAVE TO RESURRECT IT
    # if (ZOLTAN_NEEDS_ParMETIS)
    #   find_package(ParMETIS REQUIRED)
    #   if (NOT ParMETIS_FOUND)
    #     message(FATAL_ERROR  "Error: could not find ParMETIS library \
    # 	    needed by ZOLTAN\n \
    # 	   Try defining ParMETIS_DIR: \
    # 	   (-D ParMETIS_DIR:FILEPATH=<ParMETIS_install_prefix>)\n")
    #   endif (NOT ParMETIS_FOUND)
    
    #   list(APPEND ZOLTAN_LIBRARIES ${ParMETIS_LIBRARIES})
    #   message("ParMETIS_LIBRARIES --- " ${ParMETIS_LIBRARIES})
    #   message("ZOLTAN_LIBRARIES --- " ${Zoltan_LIBRARIES})
    # endif (ZOLTAN_NEEDS_ParMETIS)
    
    # if (ZOLTAN_NEEDS_PTScotch)
    #   find_package(PTScotch REQUIRED)
    #   if (PTScotch-NOTFOUND)
    #     message(FATAL_ERROR  "Error: could not find PTScotch library \
    # 	    needed by ZOLTAN\n \
    # 	   Try defining PTScotch_DIR: \
    # 	   (-D PTScotch_DIR:FILEPATH=<PTScotch_install_prefix>)\n")
    #   endif (PTScotch-NOTFOUND)
    
    #   list(APPEND ZOLTAN_LIBRARIES ${PTScotch_LIBRARIES})
    #   message("PTScotch_LIBRARIES --- " ${PTScotch_LIBRARIES})
    #   message("ZOLTAN_LIBRARIES --- " ${Zoltan_LIBRARIES})
    # endif (ZOLTAN_NEEDS_PTScotch)
    
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
  
  set(BUILD_TESTS TRUE)

  # For backward compatibility allow UnitTest_DIR for a bit more time
  find_package(UnitTest++ HINTS ${UnitTest++_DIR} ${UnitTest_DIR})
  
  if (UnitTest++_FOUND)
    set(BUILD_TESTS TRUE)
    
    # UnitTest++Config.cmake uses a weird name for storing it's include path
    # and it's not even completely right
    set(UnitTest++_INCLUDE_DIRS ${UTPP_INCLUDE_DIRS}/UnitTest++)
    
    # Also it does not set this variable
    set(UnitTest++_LIBRARIES UnitTest++)
    
  else (UnitTest++_FOUND)

    # can happen if unittest is installed using a package manager like apt
    message(STATUS "find_package did not find CMake config for UnitTest++. Find manually")
    # Find things by hand
    
    find_library(UnitTest++_LIBRARY
      NAMES UnitTest++ unittest++
      HINTS
      ${UnitTest++_DIR}
      ${UnitTest++_LIBRARY_DIR}
      ${UnitTest_DIR}
      ${UnitTest_LIBRARY_DIR}
      )
    if (NOT UnitTest++_LIBRARY)
      message(FATAL_ERROR "Error: could not find unittest++ library")
    endif(NOT UnitTest++_LIBRARY)
    get_filename_component(UnitTest++_LIBRARY_DIRS ${UnitTest++_LIBRARY} DIRECTORY)
    
    find_path(UnitTest++_INCLUDE_DIRS "UnitTest++.h" 
      HINTS 
      ${UnitTest++_DIR}
      ${UnitTest++_DIR}/src 
      ${UnitTest++_DIR}/include/unittest++
      ${UnitTest++_DIR}/include/UnitTest++
      ${UnitTest_DIR}
      ${UnitTest_DIR}/src 
      ${UnitTest_DIR}/include/unittest++
      ${UnitTest_DIR}/include/UnitTest++
      )
    if (NOT UnitTest++_INCLUDE_DIRS)
      message(FATAL_ERROR "Error: could not find UnitTest++ includes")
    endif(NOT UnitTest++_INCLUDE_DIRS)
    
    set(UnitTest++_LIBRARIES UnitTest++)

    add_library(UnitTest++ UNKNOWN IMPORTED)
    set_target_properties(UnitTest++ PROPERTIES IMPORTED_LOCATION ${UnitTest++_LIBRARY})
  endif(UnitTest++_FOUND)

  message(STATUS "### UnitTest++ library: ${UnitTest++_LIBRARIES}")
  get_target_property(UnitTest++_LOC UnitTest++ IMPORTED_LOCATION)
  if (UnitTest++_LOC)
    message(STATUS "### UnitTest++ location: ${UnitTest++_LOC}")
  else (UnitTest++_LOC)
    if (UnitTest++_LIBRARY_DIRS)
      message(STATUS "### UnitTest++ library dir: ${UnitTest++_LIBRARY_DIRS}")
    endif (UnitTest++_LIBRARY_DIRS)
  endif (UnitTest++_LOC)
  
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



