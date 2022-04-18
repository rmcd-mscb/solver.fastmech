add_test(
  NAME simple-0
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_SOURCE_DIR}/Case1.cgn.orig ${CMAKE_CURRENT_BINARY_DIR}/Case1.cgn
  )

add_test(
  NAME simple-1
  COMMAND $<TARGET_FILE:Fastmech> Case1.cgn
  )

add_test(
  NAME simple-2
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/Case1.cgn ${CMAKE_CURRENT_BINARY_DIR}/Case1.1.cgn
  )

add_test(
  NAME simple-3
  COMMAND $<TARGET_FILE:Fastmech> Case1.cgn
  )

add_test(
  NAME simple-4
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${CMAKE_CURRENT_BINARY_DIR}/Case1.cgn ${CMAKE_CURRENT_BINARY_DIR}/Case1.2.cgn
  )

# # add_test(
# #   NAME simple-5
# #   COMMAND hdf5::h5diff -c ${CMAKE_CURRENT_BINARY_DIR}/Case1.1.cgn ${CMAKE_CURRENT_BINARY_DIR}/Case1.2.cgn
# #   )


set_tests_properties(simple-1 PROPERTIES DEPENDS "simple-0")
set_tests_properties(simple-2 PROPERTIES DEPENDS "simple-1")
set_tests_properties(simple-3 PROPERTIES DEPENDS "simple-2")
set_tests_properties(simple-4 PROPERTIES DEPENDS "simple-3")