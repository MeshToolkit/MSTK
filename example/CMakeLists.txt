
include_directories(${PROJECT_SOURCE_DIR}/include)


add_executable(serial_example src/serial_example.c)
target_link_libraries(serial_example ${MSTKLIB} ${ExodusII_LIBRARIES} ${NetCDF_LIBRARIES} ${HDF5_LIBRARIES} ${METIS_LIBRARIES} ${ZOLTAN_LIBRARIES} m)
install(TARGETS serial_example DESTINATION ${PROJECT_SOURCE_DIR}/example)

if (ENABLE_PARALLEL) 
  add_executable(parallel_example src/parallel_example.c)
  target_link_libraries(parallel_example ${MSTKLIB} ${ExodusII_LIBRARIES} ${NetCDF_LIBRARIES} ${HDF5_LIBRARIES} ${METIS_LIBRARIES} ${ZOLTAN_LIBRARIES} m)
  install(TARGETS parallel_example DESTINATION ${PROJECT_SOURCE_DIR}/example)
endif ()
