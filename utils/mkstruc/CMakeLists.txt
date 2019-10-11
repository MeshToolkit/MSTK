# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


include_directories(${PROJECT_BINARY_DIR}/include)

add_executable(mkstruc src/mkstruc.c src/transform_xyz.c)
target_link_libraries(mkstruc 
                      ${MSTKLIB} 
		      ${MSTK_TPL_LIBRARIES}
		      m)

if (INSTALL_PREFIX_ARCHOS)
  install(TARGETS mkstruc DESTINATION ${INSTALL_DIR}/bin/${ARCHOS})
else ()
  install(TARGETS mkstruc DESTINATION ${INSTALL_DIR}/bin)
endif ()


