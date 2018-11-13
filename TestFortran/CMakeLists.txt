include(CTest)
include_directories(${CMAKE_Fortran_MODULE_DIRECTORY})

set(test1 GetValueTest)
add_test(NAME ${test1} COMMAND ${test1})
set(src_${test1} ${test1}.f90 testing_helpers.f90)
add_executable(${test1} ${src_${test1}})
target_link_libraries(${test1} ${bmifastmech_lib} ${EXTRA_LIBS})
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/Test1.cgn
  DESTINATION ${CMAKE_BINARY_DIR}/test
)

set(test2 conflicting_instances_test)
add_test(NAME ${test2} COMMAND ${test2})
set(src_${test2} ${test2}.f90 testing_helpers.f90)
add_executable(${test2} ${src_${test2}})
target_link_libraries(${test2} ${bmifastmech_lib} ${EXTRA_LIBS})
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/Case2.cgn
  DESTINATION ${CMAKE_BINARY_DIR}/test
)
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/Case2_1.cgn
  DESTINATION ${CMAKE_BINARY_DIR}/test
)


set(test3 irf_test)
add_test(NAME ${test3} COMMAND ${test3})
set(src_${test3} ${test3}.f90 testing_helpers.f90)
add_executable(${test3} ${src_${test3}})
target_link_libraries(${test3} ${bmifastmech_lib} ${EXTRA_LIBS})
file(
  COPY ${CMAKE_CURRENT_SOURCE_DIR}/Test1.cgn
  DESTINATION ${CMAKE_BINARY_DIR}/test
)
