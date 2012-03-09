# -*- mode: cmake -*-

#
# ZOLTAN Find Module for MSTK
# Shamelessly adapted from Amanzi open source code https://software.lanl.gov/ascem/trac
#
# Usage:
#    Control the search through ZOLTAN_DIR or setting environment variable
#    ZOLTAN_ROOT to the ZOLTAN installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    ZOLTAN_FOUND            (BOOL)       Flag indicating if ZOLTAN was found
#    ZOLTAN_INCLUDE_DIR      (PATH)       Path to the ZOLTAN include file
#    ZOLTAN_INCLUDE_DIRS     (LIST)       List of all required include files
#    ZOLTAN_LIBRARY_DIR      (PATH)       Path to the ZOLTAN library
#    ZOLTAN_LIBRARY          (FILE)       ZOLTAN library
#    ZOLTAN_LIBRARIES        (LIST)       List of all required ZOLTAN libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)

if ( ZOLTAN_LIBRARIES AND ZOLTAN_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(ZOLTAN_LIBRARIES AND ZOLTAN_INCLUDE_DIRS)

    # Cache variables
    if(ZOLTAN_DIR)
        set(ZOLTAN_DIR "${ZOLTAN_DIR}" CACHE PATH "Path to search for ZOLTAN include and library files")
    endif()

    if(ZOLTAN_INCLUDE_DIR)
        set(ZOLTAN_INCLUDE_DIR "${ZOLTAN_INCLUDE_DIR}" CACHE PATH "Path to search for ZOLTAN include files")
    endif()

    if(ZOLTAN_LIBRARY_DIR)
        set(ZOLTAN_LIBRARY_DIR "${ZOLTAN_LIBRARY_DIR}" CACHE PATH "Path to search for ZOLTAN library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) ZOLTAN_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) ZOLTAN_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(zoltan_inc_names "zoltan.h")
    if (ZOLTAN_INCLUDE_DIR)

        if (EXISTS "${ZOLTAN_INCLUDE_DIR}")

            find_path(zoltan_test_include_path
                      NAMES ${zoltan_inc_names}
                      HINTS ${ZOLTAN_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT zoltan_test_include_path)
                message(SEND_ERROR "Can not locate ${zoltan_inc_names} in ${ZOLTAN_INCLUDE_DIR}")
            endif()
            set(ZOLTAN_INCLUDE_DIR "${zoltan_test_include_path}")

        else()
            message(SEND_ERROR "ZOLTAN_INCLUDE_DIR=${ZOLTAN_INCLUDE_DIR} does not exist")
            set(ZOLTAN_INCLUDE_DIR "ZOLTAN_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(zoltan_inc_suffixes "include")
        if(ZOLTAN_DIR)

            if (EXISTS "${ZOLTAN_DIR}" )

                find_path(ZOLTAN_INCLUDE_DIR
                          NAMES ${zoltan_inc_names}
                          HINTS ${ZOLTAN_DIR}
                          PATH_SUFFIXES ${zoltan_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "ZOLTAN_DIR=${ZOLTAN_DIR} does not exist")
                 set(ZOLTAN_INCLUDE_DIR "ZOLTAN_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(ZOLTAN_INCLUDE_DIR
                      NAMES ${zoltan_inc_names}
                      PATH_SUFFIXES ${zoltan_inc_suffixes})

        endif()

    endif()

    if ( NOT ZOLTAN_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate ZOLTAN include directory")
    endif()

    # Search for libraries 
    # Search order preference:
    #  (1) ZOLTAN_LIBRARY_DIR - check existence of path AND if the library file exists
    #  (2) ZOLTAN_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(zoltan_lib_names "zoltan")
    if (ZOLTAN_LIBRARY_DIR)

        if (EXISTS "${ZOLTAN_LIBRARY_DIR}")

            find_library(ZOLTAN_LIBRARY
                         NAMES ${zoltan_lib_names}
                         HINTS ${ZOLTAN_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "ZOLTAN_LIBRARY_DIR=${ZOLTAN_LIBRARY_DIR} does not exist")
            set(ZOLTAN_LIBRARY "ZOLTAN_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND zoltan_lib_suffixes "lib" "Lib")
        if(ZOLTAN_DIR)

            if (EXISTS "${ZOLTAN_DIR}" )

                find_library(ZOLTAN_LIBRARY
                             NAMES ${zoltan_lib_names}
                             HINTS ${ZOLTAN_DIR}
                             PATH_SUFFIXES ${zoltan_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "ZOLTAN_DIR=${ZOLTAN_DIR} does not exist")
                 set(ZOLTANLIBRARY "ZOLTAN_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(ZOLTAN_LIBRARY
                         NAMES ${zoltan_lib_names}
                         PATH_SUFFIXES ${zoltan_lib_suffixes})

        endif()

    endif()

    if ( NOT ZOLTAN_LIBRARY )
        message(SEND_ERROR "Can not locate ZOLTAN library")
    endif()    

   
    # Define prerequisite packages
    set(ZOLTAN_INCLUDE_DIRS ${ZOLTAN_INCLUDE_DIR})
    set(ZOLTAN_LIBRARIES    ${ZOLTAN_LIBRARY})

   
endif(ZOLTAN_LIBRARIES AND ZOLTAN_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(ZOLTAN DEFAULT_MSG
                                  ZOLTAN_LIBRARIES
                                  ZOLTAN_INCLUDE_DIRS)

# find_package_handle_standard_args should set ZOLTAN_FOUND but it does not!
if ( ZOLTAN_LIBRARIES AND ZOLTAN_INCLUDE_DIRS)
    set(ZOLTAN_FOUND TRUE)
else()
    set(ZOLTAN_FOUND FALSE)
endif()

# Define the version

mark_as_advanced(
  ZOLTAN_INCLUDE_DIR
  ZOLTAN_INCLUDE_DIRS
  ZOLTAN_LIBRARY
  ZOLTAN_LIBRARIES
  ZOLTAN_LIBRARY_DIR
)
