
if (BUILD_TESTS)

 # Need to copy files for tests otherwise the input files aren't accessible

  if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
      copy_directory ${MSTK_SOURCE_DIR}/unittests ${MSTK_BINARY_DIR}/unittests)
  endif ()

  include_directories(${UnitTest_INCLUDE_DIRS})

  file (GLOB FILES serial/*.cc)
  foreach (FILE ${FILES})
    get_filename_component(RELFILE ${FILE} NAME) 
    if (${RELFILE} MATCHES ExodusII OR $(RELFILE) MATCHES NemesisI)
      if (ENABLE_ExodusII)
        set (serialtestfiles ${serialtestfiles} ${FILE})
      endif ()
    else ()
      set (serialtestfiles ${serialtestfiles} ${FILE})
    endif ()
  endforeach ()

  add_executable(mstk_serial_unit_tests ${serialtestfiles})

  target_link_libraries(mstk_serial_unit_tests 
    ${MSTKLIB} 
    ${METIS_LIBRARIES} 
    ${ZOLTAN_LIBRARIES} 
    ${ExodusII_LIBRARIES} 
    ${NetCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${UNITTEST_LIBRARY})

  add_test(NAME serial_unittests COMMAND mstk_serial_unit_tests)



  if (ENABLE_PARALLEL)

    file (GLOB FILES parallel/4proc/*.cc)
    foreach (FILE ${FILES})
      get_filename_component(RELFILE ${FILE} NAME) 
      if (${RELFILE} MATCHES ExodusII OR $(RELFILE) MATCHES NemesisI)
        if (ENABLE_ExodusII)
          set (par4testfiles ${par4testfiles} ${FILE})
        endif ()
      else ()
        set (par4testfiles ${par4testfiles} ${FILE})
      endif ()
    endforeach ()

    add_executable(mstk_parallel_4proc_unit_tests ${par4testfiles})

    target_link_libraries(mstk_parallel_4proc_unit_tests 
      ${MSTKLIB} 
      ${METIS_LIBRARIES} 
      ${ZOLTAN_LIBRARIES} 
      ${ExodusII_LIBRARIES} 
      ${NetCDF_LIBRARIES}
      ${HDF5_LIBRARIES}
      ${UNITTEST_LIBRARY})

    add_test(NAME parallel_unittests_4proc 
      COMMAND ${MPIEXEC} -n 4 mstk_parallel_4proc_unit_tests)
    


    file (GLOB FILES parallel/8proc/*.cc)
    foreach (FILE ${FILES})
      get_filename_component(RELFILE ${FILE} NAME) 
      if (${RELFILE} MATCHES ExodusII OR $(RELFILE) MATCHES NemesisI)
        if (ENABLE_ExodusII)
          set (par8testfiles ${par8testfiles} ${FILE})
        endif ()
      else ()
        set (par8testfiles ${par8testfiles} ${FILE})
      endif ()
    endforeach ()

    add_executable(mstk_parallel_8proc_unit_tests ${par8testfiles})

    target_link_libraries(mstk_parallel_8proc_unit_tests 
      ${MSTKLIB} 
      ${METIS_LIBRARIES} 
      ${ZOLTAN_LIBRARIES} 
      ${ExodusII_LIBRARIES} 
      ${NetCDF_LIBRARIES}
      ${HDF5_LIBRARIES}
      ${UNITTEST_LIBRARY})

    add_test(NAME parallel_unittests_8proc 
      COMMAND ${MPIEXEC} -n 8 mstk_parallel_8proc_unit_tests)

    endif ()

endif ()
