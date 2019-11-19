# Copyright 2019 Triad National Security, LLC. All rights reserved.
# 
# This file is part of the MSTK project. Please see the license file at
# the root of this repository or at
# https://github.com/MeshToolkit/MSTK/blob/master/LICENSE


if (BUILD_TESTS)

  # We need to do some tests with multi-threading to see if the code is thread-safe
  find_package(OpenMP)
  if (OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif(OPENMP_FOUND)
  
  # Need to copy files for tests otherwise the input files aren't accessible

  if (NOT (${PROJECT_SOURCE_DIR} EQUAL ${PROJECT_BINARY_DIR}))
    execute_process(COMMAND ${CMAKE_COMMAND} -E
      copy_directory ${PROJECT_SOURCE_DIR}/unittests ${PROJECT_BINARY_DIR}/unittests)
  endif ()

  file (GLOB FILES serial/*.cc)
  foreach (FILE ${FILES})
    get_filename_component(RELFILE ${FILE} NAME) 
    if (${RELFILE} MATCHES ExodusII OR $(RELFILE) MATCHES NemesisI)
      
      if (ENABLE_ExodusII)
        set (serialtestfiles ${serialtestfiles} ${FILE})
      endif (ENABLE_ExodusII)

    elseif (${RELFILE} EQUAL "test_Markers.cc")

      if (OPENMP_FOUND)
	set (serialtestfiles ${serialtestfiles} ${FILE})
      endif (OPENMP_FOUND)

    else (${RELFILE} MATCHES ExodusII OR $(RELFILE) MATCHES NemesisI)

      set (serialtestfiles ${serialtestfiles} ${FILE})

    endif (${RELFILE} MATCHES ExodusII OR $(RELFILE) MATCHES NemesisI)
  endforeach ()

  add_executable(mstk_serial_unit_tests ${serialtestfiles})

  target_include_directories(mstk_serial_unit_tests
    PRIVATE ${UnitTest++_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include)
  target_link_libraries(mstk_serial_unit_tests 
    ${MSTKLIB} 
    ${METIS_LIBRARIES} 
    ${ZOLTAN_LIBRARIES} 
    ${ExodusII_LIBRARIES} 
    ${NetCDF_LIBRARIES}
    ${HDF5_LIBRARIES}
    ${UnitTest++_LIBRARIES})

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

    target_include_directories(mstk_parallel_4proc_unit_tests
      PRIVATE ${UnitTest++_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include)
    target_link_libraries(mstk_parallel_4proc_unit_tests 
      ${MSTKLIB} 
      ${METIS_LIBRARIES} 
      ${ZOLTAN_LIBRARIES} 
      ${ExodusII_LIBRARIES} 
      ${NetCDF_LIBRARIES}
      ${HDF5_LIBRARIES}
      ${UnitTest++_LIBRARIES})

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

    target_include_directories(mstk_parallel_8proc_unit_tests
      PRIVATE ${UnitTest++_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include)
    target_link_libraries(mstk_parallel_8proc_unit_tests 
      ${MSTKLIB} 
      ${METIS_LIBRARIES} 
      ${ZOLTAN_LIBRARIES} 
      ${ExodusII_LIBRARIES} 
      ${NetCDF_LIBRARIES}
      ${HDF5_LIBRARIES}
      ${UnitTest++_LIBRARIES})

    add_test(NAME parallel_unittests_8proc 
      COMMAND ${MPIEXEC} -n 8 mstk_parallel_8proc_unit_tests)

    endif ()

endif ()
