
include_directories(${PROJECT_SOURCE_DIR}/include)


add_executable(chkmesh src/chkmesh.c src/geom.c)
target_link_libraries(chkmesh ${MSTKLIB} ${EXODUS_LIBS_ALL} ${NETCDF_LIBS_ALL})

# Would be nice to specify it relative to the path of this CMakeLists.txt file

install(TARGETS chkmesh DESTINATION ${PROJECT_SOURCE_DIR}/utils/chkmesh/bin/${ARCHOS})

