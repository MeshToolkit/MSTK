
include_directories(${PROJECT_SOURCE_DIR}/include)


add_executable(meshconvert src/meshconvert.c)
target_link_libraries(meshconvert ${MSTKLIB} ${ExodusII_LIBRARIES} ${NetCDF_LIBRARIES} ${METIS_LIBRARIES} ${ZOLTAN_LIBRARIES} m)

# Would be nice to specify it relative to the path of this CMakeLists.txt file

install(TARGETS meshconvert DESTINATION ${PROJECT_SOURCE_DIR}/utils/meshconvert/bin/${ARCHOS})

