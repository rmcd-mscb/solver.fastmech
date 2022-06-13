project(test1 Fortran)

configure_file(Case1.cgn.orig Case1.cgn.orig COPYONLY)

find_program(PWSH NAMES pwsh powershell)

add_executable(test1
  test1.f90
  ../../src/UnlinkHDF5Mod.f90
)

target_link_libraries(test1 PRIVATE hdf5::hdf5_fortran-static)

# Don't think this should be necessary
set_target_properties (hdf5::hdf5_fortran-static
  PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES ${HDF5_INCLUDE_DIRS}/static
)

# POST_BUILD
if (WIN32)
  # copy reqd dlls to test1 directory
  add_custom_command(TARGET test1 POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:hdf5::hdf5-shared>"  "$<TARGET_FILE_DIR:test1>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:zlib-shared>"        "$<TARGET_FILE_DIR:test1>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$<TARGET_FILE:szip-shared>"        "$<TARGET_FILE_DIR:test1>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$ENV{IFORT_COMPILER21}\\redist\\intel64_win\\compiler\\$<IF:$<CONFIG:Release>,libifcoremd.dll,libifcoremdd.dll>" "$<TARGET_FILE_DIR:test1>"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different "$ENV{IFORT_COMPILER21}\\redist\\intel64_win\\compiler\\$<IF:$<CONFIG:Release>,libmmd.dll,libmmdd.dll>"           "$<TARGET_FILE_DIR:test1>"
    )
endif()

# test
add_test(test.TestCleanHDF5.1 test1)

if (PWSH)
  # THIS DOESN'T WORK
  # set(groups iRIC3D BaseIterativeData)
  # foreach(group ${groups})
  #   add_test(NAME test.TestCleanHDF5.${group}.FAIL
  #     COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String Group ${group} ) -ne $null ) { exit 1 } ]=]
  #     COMMAND_EXPAND_LISTS
  #   )
  #   set_tests_properties(test.TestCleanHDF5.${group}.FAIL PROPERTIES WILL_FAIL TRUE) 
  # endforeach()

  # expected to fail ( Case1.cgn.orig)
  add_test(NAME test.unlink_solutions.BaseIterativeData.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "BaseIterativeData"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.BaseIterativeData.FAIL PROPERTIES WILL_FAIL TRUE)

  add_test(NAME test.unlink_solutions.FlowCellSolution1.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowCellSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowCellSolution1.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.FlowIFaceSolution1.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowIFaceSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowIFaceSolution1.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.FlowJFaceSolution1.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowJFaceSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowJFaceSolution1.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.FlowSolution1.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowSolution1.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.GridCoordinatesForSolution1.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "GridCoordinatesForSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.GridCoordinatesForSolution1.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.ParticleGroupSolution1.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "ParticleGroupSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.ParticleGroupSolution1.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.PolydataSolution1.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "PolydataSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.PolydataSolution1.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.FlowCellSolutionPointers.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowCellSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowCellSolutionPointers.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.FlowIFaceSolutionPointers.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowIFaceSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowIFaceSolutionPointers.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.FlowJFaceSolutionPointers.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowJFaceSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowJFaceSolutionPointers.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.FlowSolutionPointers.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "FlowSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.FlowSolutionPointers.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.GridCoordinatesPointers.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "GridCoordinatesPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.GridCoordinatesPointers.FAIL PROPERTIES WILL_FAIL TRUE) 

  add_test(NAME test.unlink_solutions.iRIC3D.FAIL
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn.orig | Select-String 'Group "iRIC3D"' ) -ne $null ) { exit 1 } ]=]
           )
  set_tests_properties(test.unlink_solutions.iRIC3D.FAIL PROPERTIES WILL_FAIL TRUE) 


  # expected to succeed (Case1.cgn)
  add_test(NAME test.unlink_solutions.BaseIterativeData
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "BaseIterativeData"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowCellSolution1
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowCellSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowIFaceSolution1
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowIFaceSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowJFaceSolution1
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowJFaceSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowSolution1
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.GridCoordinatesForSolution1
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "GridCoordinatesForSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.ParticleGroupSolution1
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "ParticleGroupSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.PolydataSolution1
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "PolydataSolution1"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowCellSolutionPointers
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowCellSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowIFaceSolutionPointers
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowIFaceSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowJFaceSolutionPointers
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowJFaceSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.FlowSolutionPointers
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "FlowSolutionPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.GridCoordinatesPointers
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "GridCoordinatesPointers"' ) -ne $null ) { exit 1 } ]=]
           )
  add_test(NAME test.unlink_solutions.iRIC3D
           COMMAND ${PWSH} -c [=[ if ( ( $<TARGET_FILE_DIR:hdf5::h5diff>/h5dump -H Case1.cgn | Select-String 'Group "iRIC3D"' ) -ne $null ) { exit 1 } ]=]
           )
endif()
