
if (BUILD_TESTS)

 # Need to copy files for tests otherwise the input files aren't accessible

  if (NOT (${MSTK_SOURCE_DIR} EQUAL ${MSTK_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
      copy_directory ${MSTK_SOURCE_DIR}/unittests ${MSTK_BINARY_DIR}/unittests)
  endif ()

  include_directories(${UnitTest_INCLUDE_DIRS})

  file (GLOB FILES serial/*.cc)
  add_executable(mstk_serial_unit_tests ${FILES})
  target_link_libraries(mstk_serial_unit_tests ${MSTKLIB} ${UNITTEST_LIBRARY} m)
  add_test(NAME serial_unittests COMMAND mstk_serial_unit_tests)

  if (ENABLE_PARALLEL) 
    file (GLOB FILES parallel/4proc/*.cc)
    add_executable(mstk_parallel_4proc_unit_tests ${FILES})
    target_link_libraries(mstk_parallel_4proc_unit_tests ${MSTKLIB} ${METIS_LIBRARY} ${UNITTEST_LIBRARY})
  add_test(NAME parallel_unittests COMMAND mpirun -n 4 mstk_parallel_4proc_unit_tests)
  endif ()

endif ()
