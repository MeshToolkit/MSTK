
include_directories(${PROJECT_SOURCE_DIR}/include)


add_executable(vizpart src/vizpart.c)
target_link_libraries(vizpart 
                      ${MSTKLIB} 
                      ${MSTK_TPL_LIBRARIES} 
                      m)

install(TARGETS vizpart DESTINATION ${PROJECT_SOURCE_DIR}/utils/vizpart/bin/${ARCHOS})

