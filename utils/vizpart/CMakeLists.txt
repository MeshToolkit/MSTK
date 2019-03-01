
include_directories(${PROJECT_BINARY_DIR}/include)


add_executable(vizpart src/vizpart.c)
target_link_libraries(vizpart 
                      ${MSTKLIB} 
                      ${MSTK_TPL_LIBRARIES} 
                      m)

if (INSTALL_DIR) 
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS vizpart DESTINATION ${INSTALL_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS vizpart DESTINATION ${INSTALL_DIR}/bin)
  endif ()
else ()
  if (INSTALL_PREFIX_ARCHOS)
    install(TARGETS vizpart DESTINATION ${PROJECT_SOURCE_DIR}/bin/${ARCHOS})
  else ()
    install(TARGETS vizpart DESTINATION ${PROJECT_SOURCE_DIR}/bin)
  endif ()
endif ()

