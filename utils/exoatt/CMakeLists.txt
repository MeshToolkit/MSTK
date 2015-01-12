
include_directories(${PROJECT_SOURCE_DIR}/include)

add_executable(exoatt src/exoatt.c)
target_link_libraries(exoatt 
                      ${MSTKLIB} 
		      ${MSTK_TPL_LIBRARIES}
		      m)

# Would be nice to specify it relative to the path of this CMakeLists.txt file

if (INSTALL_DIR) 
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS exoatt DESTINATION ${INSTALL_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS exoatt DESTINATION ${INSTALL_DIR}/bin)
  endif ()
else ()
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS exoatt DESTINATION ${PROJECT_SOURCE_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS exoatt DESTINATION ${PROJECT_SOURCE_DIR}/bin)
  endif ()
endif ()


