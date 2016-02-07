
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(mkstruc src/mkstruc.c src/transform_xyz.c)
target_link_libraries(mkstruc 
                      ${MSTKLIB} 
		      ${MSTK_TPL_LIBRARIES}
		      m)

# Would be nice to specify it relative to the path of this CMakeLists.txt file

if (INSTALL_DIR) 
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS mkstruc DESTINATION ${INSTALL_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS mkstruc DESTINATION ${INSTALL_DIR}/bin)
  endif ()
else ()
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS mkstruc DESTINATION ${PROJECT_SOURCE_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS mkstruc DESTINATION ${PROJECT_SOURCE_DIR}/bin)
  endif ()
endif ()


