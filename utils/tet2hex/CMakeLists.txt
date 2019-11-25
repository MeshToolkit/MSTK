# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


add_executable(tet2hex src/tet2hex.c)
target_link_libraries(tet2hex ${MSTKLIB})
target_include_directories(tet2hex PRIVATE ${PROJECT_BINARY_DIR}/include)


install(TARGETS tet2hex DESTINATION bin)


