# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE



add_executable(exoatt src/exoatt.c)
target_link_libraries(exoatt ${MSTKLIB})
target_include_directories(exoatt PRIVATE ${PROJECT_BINARY_DIR}/include)

install(TARGETS exoatt DESTINATION bin)

# Unittests

if (ENABLE_Tests)
  add_subdirectory(unittests)
endif()
