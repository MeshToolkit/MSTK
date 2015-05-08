
if (BUILD_TESTS)

 # Need to copy files for tests otherwise the input files aren't accessible

  if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
      copy_directory ${MSTK_SOURCE_DIR}/utils/exoatt/test ${MSTK_BINARY_DIR}/utils/exoatt/test)
  endif ()

  include_directories(${UnitTest_INCLUDE_DIRS})

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

  target_link_libraries(exoatt_unit_tests 
    ${MSTKLIB} 
    ${METIS_LIBRARIES} 
    ${ZOLTAN_LIBRARIES} 
    ${ExodusII_LIBRARIES} 
    ${NetCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${UNITTEST_LIBRARY})

  add_test(NAME exoatt_unittests COMMAND exoatt_unit_tests)


endif(BUILD_TESTS)
