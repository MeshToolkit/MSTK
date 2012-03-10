
include_directories(${PROJECT_SOURCE_DIR}/include)


add_executable(serial_example src/serial_example.c)
target_link_libraries(serial_example ${MSTKLIB} m)
install(TARGETS serial_example DESTINATION ${PROJECT_SOURCE_DIR}/example)
add_executable(test_metis src/test_metis.c)
target_link_libraries(test_metis ${MSTKLIB} m)
install(TARGETS test_metis DESTINATION ${PROJECT_SOURCE_DIR}/example)


if (ENABLE_PARALLEL)
  add_executable(parallel_example src/parallel_example.c)
  target_link_libraries(parallel_example ${MSTKLIB} ${METIS_LIBRARY} m)
  install(TARGETS parallel_example DESTINATION ${PROJECT_SOURCE_DIR}/example)

  add_executable(test_zoltan src/test_zoltan.c)
  target_link_libraries(test_zoltan ${MSTKLIB} ${Zoltan_LIBRARY} m)
  install(TARGETS test_zoltan DESTINATION ${PROJECT_SOURCE_DIR}/example)

endif ()
