# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE



add_executable(exoatt src/exoatt.c)

target_include_directories(exoatt PRIVATE ${PROJECT_BINARY_DIR}/include)
target_link_libraries(exoatt 
                      ${MSTKLIB} 
		      ${MSTK_TPL_LIBRARIES}
		      m)

if (INSTALL_PREFIX_ARCHOS)
  install(TARGETS exoatt DESTINATION ${INSTALL_DIR}/bin/${ARCHOS})
else ()
  install(TARGETS exoatt DESTINATION ${INSTALL_DIR}/bin)
endif ()

# Unittests

if (ENABLE_Tests)
  add_subdirectory(unittests)
endif()
