# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE

project (OPTMESH)

cmake_minimum_required(VERSION 2.8)

set(ARCHOS ${CMAKE_SYSTEM_PROCESSOR}_${CMAKE_SYSTEM_NAME})

#
# MSTK Libraries and includes
#

if (NOT MSTK_DIR)
  message(FATAL_ERROR "Error: MSTK top level installation dir must be defined")
endif()

# Look for the cmake configuration of the MSTK build (if it exists)

find_file(MSTK_CONFIG_FILE MSTKConfig.cmake PATH ${MSTK_DIR}/share/cmake)


if (MSTK_CONFIG_FILE) 

  # if we found the MSTK config file get all the required info from it

  include(${MSTK_CONFIG_FILE})

  find_library(MSTK_LIBRARY ${MSTK_LIBRARY_NAME} PATH ${MSTK_LIBRARY_DIR})
  if (NOT MSTK_LIBRARY)
    message(FATAL_ERROR "Could not find MSTK_LIBRARY with name ${MSTK_LIBRARY_NAME}")
  endif()

else ()

  # Did not find the MSTK cmake configuration file
  # We have to discover all the libraries and includes ourselves
  #
  # To keep example simple, just quit

  message(FATAL_ERROR "Could not find MSTK configuration file")

endif()

# Now gather the source files for this application, 
# and link executable

file (GLOB appsrcfiles src/*.c)

add_executable(myapp ${appsrcfiles})
include_directories(${PROJECT_SOURCE_DIR}/include ${MSTK_INCLUDE_DIR})
target_link_libraries(myapp ${MSTK_LIBRARY} ${MSTK_TPL_LIBRARIES} z m)


# install

if (CMAKE_INSTALL_PREFIX)
  install(TARGETS myapp DESTINATION ${CMAKE_INSTALL_PREFIX}/bin/${ARCHOS}/myapp)
else ()
  install(TARGETS myapp DESTINATION ${PROJECT_SOURCE_DIR}/bin/${ARCHOS}/myapp)
endif ()


