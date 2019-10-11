# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


include_directories(${PROJECT_BINARY_DIR}/include)

add_executable(tri2quad src/tri2quad.c)
target_link_libraries(tri2quad
                      ${MSTKLIB} 
		      ${MSTK_TPL_LIBRARIES}
		      m)

if (INSTALL_PREFIX_ARCHOS)
  install(TARGETS tri2quad DESTINATION ${INSTALL_DIR}/bin/${ARCHOS})
else ()
  install(TARGETS tri2quad DESTINATION ${INSTALL_DIR}/bin)
endif ()


