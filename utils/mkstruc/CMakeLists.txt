# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


add_executable(mkstruc src/mkstruc.c src/transform_xyz.c)
target_link_libraries(mkstruc ${MSTKLIB})
target_include_directories(mkstruc PRIVATE ${PROJECT_BINARY_DIR}/include)


install(TARGETS mkstruc DESTINATION bin)



