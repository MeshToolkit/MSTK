
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(tri2quad src/tri2quad.c)
target_link_libraries(tri2quad
                      ${MSTKLIB} 
		      ${MSTK_TPL_LIBRARIES}
		      m)

# Would be nice to specify it relative to the path of this CMakeLists.txt file

if (INSTALL_DIR) 
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS tri2quad DESTINATION ${INSTALL_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS tri2quad DESTINATION ${INSTALL_DIR}/bin)
  endif ()
else ()
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS tri2quad DESTINATION ${PROJECT_SOURCE_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS tri2quad DESTINATION ${PROJECT_SOURCE_DIR}/bin)
  endif ()
endif ()


