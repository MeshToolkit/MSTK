# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


add_executable(vizpart src/vizpart.c)
target_link_libraries(vizpart ${MSTKLIB})
target_include_directories(vizpart PRIVATE ${PROJECT_BINARY_DIR}/include)

install(TARGETS vizpart DESTINATION bin)


