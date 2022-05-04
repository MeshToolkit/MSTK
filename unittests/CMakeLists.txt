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

  set(serialtestfiles
    serial/Main.cc
    serial/test_MR_Edges_Hex_FN.cc
    serial/test_Attributes.cc
    serial/test_MR_Edges_Polyhed_FN.cc
    serial/test_MR_Edges_Prism_FN.cc
    serial/test_MR_Edges_Pyramid_FN.cc
    serial/test_MR_Edges_Tet_FN.cc
    serial/test_ME_Split.cc
    serial/test_MR_Rem_Face_FNR3R4.cc
    serial/test_ME_Split_SimplexMesh.cc
    serial/test_MR_Vertices_Hex_FN.cc
    serial/test_MF_Rem_Edge_F1F3.cc
    serial/test_MR_Vertices_Polyhed_FN.cc
    serial/test_MFs_Dihedral_Angle.cc
    serial/test_MR_Vertices_Prism_FN.cc
    serial/test_MF_Split.cc
    serial/test_MR_Vertices_Pyramid_FN.cc
    serial/test_MF_Split_SimplexMesh.cc
    serial/test_MR_Vertices_Tet_FN.cc
    serial/test_MF_Split_with_Edge.cc
    )
  
  if (ENABLE_ExodusII)
    list(APPEND serialtestfiles
      serial/test_Read_Write_ExodusII.cc
      serial/test_ME_Collapse.cc
      serial/test_Graph_Read_ExodusII.cc
      )

    if (OPENMP_FOUND)
      list(APPEND serialtestfiles
	serial/test_Markers.cc)
    endif (OPENMP_FOUND)
  endif (ENABLE_ExodusII)

  add_executable(mstk_serial_unit_tests ${serialtestfiles})

  target_include_directories(mstk_serial_unit_tests
    PRIVATE ${UnitTest++_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include)
  target_link_libraries(mstk_serial_unit_tests PRIVATE
    ${MSTKLIB}  ${UnitTest++_LIBRARIES})

  add_test(NAME serial_unittests COMMAND mstk_serial_unit_tests)



  if (ENABLE_PARALLEL)

    set(par4testfiles
      parallel/4proc/Main.cc
      parallel/4proc/Test_partition.cc
      parallel/4proc/Test_RenumberGlobalIDs.cc
      parallel/4proc/Test_UpdAtt.cc
      parallel/4proc/Test_UpdVtxCoords.cc
      parallel/4proc/Test_weave2D.cc
      )

    if (ENABLE_ExodusII)
      list(APPEND par4testfiles
	parallel/4proc/test_ME_Collapse_4proc.cc
	parallel/4proc/Test_partition2D_ExodusII.cc)
    endif (ENABLE_ExodusII)

    add_executable(mstk_parallel_4proc_unit_tests ${par4testfiles})

    target_include_directories(mstk_parallel_4proc_unit_tests
      PRIVATE ${UnitTest++_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include)
    target_link_libraries(mstk_parallel_4proc_unit_tests PRIVATE
      ${MSTKLIB} ${UnitTest++_LIBRARIES})

    add_test(NAME parallel_unittests_4proc 
      COMMAND ${MPIEXEC} -n 4 mstk_parallel_4proc_unit_tests)
    


    set(par8testfiles
      parallel/8proc/Main.cc
      parallel/8proc/Test_partition3D.cc
      parallel/8proc/Test_weave3D.cc
      )
    if (ENABLE_ExodusII)
      list(APPEND par8testfiles
	parallel/8proc/Test_partition3D_ExodusII.cc)
    endif ()

    add_executable(mstk_parallel_8proc_unit_tests ${par8testfiles})

    target_include_directories(mstk_parallel_8proc_unit_tests
      PRIVATE ${UnitTest++_INCLUDE_DIRS} ${PROJECT_BINARY_DIR}/include)
    target_link_libraries(mstk_parallel_8proc_unit_tests PRIVATE
      ${MSTKLIB} ${UnitTest++_LIBRARIES})

    add_test(NAME parallel_unittests_8proc 
      COMMAND ${MPIEXEC} -n 8 mstk_parallel_8proc_unit_tests)

    endif ()

endif ()
