# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


if (BUILD_TESTS)

 # Need to copy files for tests otherwise the input files aren't accessible

  if (NOT (${PROJECT_SOURCE_DIR} EQUAL ${PROJECT_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
      copy_directory ${PROJECT_SOURCE_DIR}/utils/exoatt/unittests ${PROJECT_BINARY_DIR}/utils/exoatt/unittests)
  endif ()

  file (GLOB FILES *.cc)
  foreach (FILE ${FILES})
    get_filename_component(RELFILE ${FILE} NAME) 
    if (${RELFILE} MATCHES ExodusII OR $(RELFILE) MATCHES NemesisI)
      if (ENABLE_ExodusII)
        set (testfiles ${testfiles} ${FILE})
      endif ()
    else ()
      set (testfiles ${testfiles} ${FILE})
    endif ()
  endforeach ()

  add_executable(exoatt_unit_tests ${testfiles})

  target_include_directories(exoatt_unit_tests
    PRIVATE ${UnitTest++_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include)

  target_link_libraries(exoatt_unit_tests 
    ${MSTKLIB} 
    ${METIS_LIBRARIES} 
    ${ZOLTAN_LIBRARIES} 
    ${ExodusII_LIBRARIES} 
    ${NetCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    UnitTest++)

  add_test(NAME exoatt_unittests COMMAND exoatt_unit_tests)


endif(BUILD_TESTS)
