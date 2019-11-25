# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


add_executable(meshconvert src/meshconvert.c)
target_link_libraries(meshconvert ${MSTKLIB})
target_include_directories(meshconvert PRIVATE ${PROJECT_BINARY_DIR}/include)

install(TARGETS meshconvert DESTINATION bin)


