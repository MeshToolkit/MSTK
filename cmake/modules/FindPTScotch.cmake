# -*- mode: cmake -*-

#
# PTScotch Find Module for MSTK
#
# Very loosely adapted from FindPTScotch modules found at:
# https://github.com/OPM/opm-common/blob/master/cmake/Modules/FindPTScotch.cmake
# https://github.com/dune-project/dune-common/blob/master/cmake/modules/FindPTScotch.cmake
#
# Usage:
#    Control the search through PTScotch_DIR or setting environment variable
#    PTScotch_ROOT to the PTScotch installation prefix.
#
#    This module does not search default paths! 
#
#    Following variables are set:
#    PTScotch_FOUND            (BOOL)     Flag indicating if PTScotch was found
#    PTScotch_INCLUDE_DIR      (PATH)     Path to the PTScotch include file
#    PTScotch_INCLUDE_DIRS     (LIST)     List of all required include files
#    PTScotch_LIBRARY_DIR      (PATH)     Path to the PTScotch library
#    PTScotch_LIBRARY          (FILE)     PTScotch library
#    PTScotch_LIBRARIES        (LIST)     List of all required PTScotch libraries
#
# #############################################################################

# Standard CMake modules see CMAKE_ROOT/Modules
include(FindPackageHandleStandardArgs)

# Amanzi CMake functions see <root>/tools/cmake for source
include(PrintVariable)

if ( PTScotch_LIBRARIES AND PTScotch_INCLUDE_DIRS )

    # Do nothing. Variables are set. No need to search again

else(PTScotch_LIBRARIES AND PTScotch_INCLUDE_DIRS)

    # Cache variables
    if(PTScotch_DIR)
        set(PTScotch_DIR "${PTScotch_DIR}" CACHE PATH "Path to search for PTScotch include and library files")
    endif()

    if(PTScotch_INCLUDE_DIR)
        set(PTScotch_INCLUDE_DIR "${PTScotch_INCLUDE_DIR}" CACHE PATH "Path to search for PTScotch include files")
    endif()

    if(PTScotch_LIBRARY_DIR)
        set(PTScotch_LIBRARY_DIR "${PTScotch_LIBRARY_DIR}" CACHE PATH "Path to search for PTScotch library files")
    endif()

    
    # Search for include files
    # Search order preference:
    #  (1) PTScotch_INCLUDE_DIR - check existence of path AND if the include files exist
    #  (2) PTScotch_DIR/<include>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(ptscotch_inc_names "ptscotch.h")
    if (PTScotch_INCLUDE_DIR)

        if (EXISTS "${PTScotch_INCLUDE_DIR}")

            find_path(ptscotch_test_include_path
                      NAMES ${ptscotch_inc_names}
                      HINTS ${PTScotch_INCLUDE_DIR}
                      NO_DEFAULT_PATH)
            if(NOT ptscotch_test_include_path)
                message(SEND_ERROR "Can not locate ${ptscotch_inc_names} in ${PTScotch_INCLUDE_DIR}")
            endif()
            set(PTScotch_INCLUDE_DIR "${ptscotch_test_include_path}")

        else()
            message(SEND_ERROR "PTScotch_INCLUDE_DIR=${PTScotch_INCLUDE_DIR} does not exist")
            set(PTScotch_INCLUDE_DIR "PTScotch_INCLUDE_DIR-NOTFOUND")
        endif()

   else() 

        set(ptscotch_inc_suffixes "include")
        if(PTScotch_DIR)

            if (EXISTS "${PTScotch_DIR}" )

                find_path(PTScotch_INCLUDE_DIR
                          NAMES ${ptscotch_inc_names}
                          HINTS ${PTScotch_DIR}
                          PATH_SUFFIXES ${ptscotch_inc_suffixes}
                          NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "PTScotch_DIR=${PTScotch_DIR} does not exist")
                 set(PTScotch_INCLUDE_DIR "PTScotch_INCLUDE_DIR-NOTFOUND")
            endif()    


        else()

            find_path(PTScotch_INCLUDE_DIR
                      NAMES ${ptscotch_inc_names}
                      PATH_SUFFIXES ${ptscotch_inc_suffixes})

        endif()

    endif()

    if ( NOT PTScotch_INCLUDE_DIR )
        message(SEND_ERROR "Can not locate PTScotch include directory")
    endif()

    # Search for libraries scotch, ptscotch and ptscotcherr
    # Search order preference:
    #  (1) PTScotch_LIBRARY_DIR - check existence of path AND if the library file exists
    #  (2) PTScotch_DIR/<lib,Lib>
    #  (3) Default CMake paths See cmake --html-help=out.html file for more information.
    #
    set(scotch_lib_names "scotch")
    set(ptscotch_lib_names "ptscotch")
    set(ptscotcherr_lib_names "ptscotcherr")

    if (PTScotch_LIBRARY_DIR)

        if (EXISTS "${PTScotch_LIBRARY_DIR}")

            find_library(Scotch_LIBRARY
                         NAMES ${scotch_lib_names}
                         HINTS ${PTScotch_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
            find_library(PTScotch_LIBRARY
                         NAMES ${ptscotch_lib_names}
                         HINTS ${PTScotch_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
            find_library(PTScotchErr_LIBRARY
                         NAMES ${ptscotcherr_lib_names}
                         HINTS ${PTScotch_LIBRARY_DIR}
                         NO_DEFAULT_PATH)
        else()
            message(SEND_ERROR "PTScotch_LIBRARY_DIR=${PTScotch_LIBRARY_DIR} does not exist")
            set(PTScotch_LIBRARY "PTScotch_LIBRARY-NOTFOUND")
        endif()

    else() 

        list(APPEND ptscotch_lib_suffixes "lib")
        if (PTScotch_DIR)

            if (EXISTS "${PTScotch_DIR}" )

                find_library(Scotch_LIBRARY
                             NAMES ${scotch_lib_names}
                             HINTS ${PTScotch_DIR}
                             PATH_SUFFIXES ${ptscotch_lib_suffixes}
                             NO_DEFAULT_PATH)
                find_library(PTScotch_LIBRARY
                             NAMES ${ptscotch_lib_names}
                             HINTS ${PTScotch_DIR}
                             PATH_SUFFIXES ${ptscotch_lib_suffixes}
                             NO_DEFAULT_PATH)
                find_library(PTScotchErr_LIBRARY
                             NAMES ${ptscotcherr_lib_names}
                             HINTS ${PTScotch_DIR}
                             PATH_SUFFIXES ${ptscotch_lib_suffixes}
                             NO_DEFAULT_PATH)

            else()
                 message(SEND_ERROR "PTScotch_DIR=${PTScotch_DIR} does not exist")
                 set(PTScotchLIBRARY "PTScotch_LIBRARY-NOTFOUND")
            endif()    


        else()

            find_library(Scotch_LIBRARY
                         NAMES ${scotch_lib_names}
                         PATH_SUFFIXES ${ptscotch_lib_suffixes})
            find_library(PTScotch_LIBRARY
                         NAMES ${ptscotch_lib_names}
                         PATH_SUFFIXES ${ptscotch_lib_suffixes})
            find_library(PTScotchErr_LIBRARY
                         NAMES ${ptscotcherr_lib_names}
                         PATH_SUFFIXES ${ptscotch_lib_suffixes})

        endif()

    endif()

    if ( NOT Scotch_LIBRARY )
      message(SEND_ERROR "Cannot locate Scotch library - libscotch")
    endif()
    if ( NOT PTScotch_LIBRARY )
        message(SEND_ERROR "Cannot locate PTScotch library - libptscotch")
    endif()    
    if ( NOT PTScotchErr_LIBRARY )
        message(SEND_ERROR "Cannot locate PTScotchErr library - libptscotcherr")
    endif()    
    

   
    # Define prerequisite packages
    set(PTScotch_INCLUDE_DIRS ${PTScotch_INCLUDE_DIR})
    set(PTScotch_LIBRARIES    ${Scotch_LIBRARY} ${PTScotch_LIBRARY} ${PTScotchErr_LIBRARY})

   
endif(PTScotch_LIBRARIES AND PTScotch_INCLUDE_DIRS )    

# Send useful message if everything is found
find_package_handle_standard_args(PTScotch DEFAULT_MSG
                                  PTScotch_LIBRARIES
                                  PTScotch_INCLUDE_DIRS)

# find_package_handle_standard_args should set PTScotch_FOUND but it does not!
if ( PTScotch_LIBRARIES AND PTScotch_INCLUDE_DIRS)
    set(PTScotch_FOUND TRUE)
else()
    set(PTScotch_FOUND FALSE)
endif()

# Define the version

mark_as_advanced(
  PTScotch_INCLUDE_DIR
  PTScotch_INCLUDE_DIRS
  PTScotch_LIBRARY
  PTScotch_LIBRARIES
  PTScotch_LIBRARY_DIR
)
