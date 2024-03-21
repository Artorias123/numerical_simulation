function(cc_library)
  set(options "")
  set(oneValueArgs NAME)
  set(multiValueArgs SRCS DEPS HDRS)
  cmake_parse_arguments(CC_LIB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
#  message("${CC_LIB_NAME}, ${CC_LIB_SRCS}, ${CC_LIB_DEPS}, ${CC_LIB_HDRS}")
  if(CC_LIB_SRCS)
    add_library(${CC_LIB_NAME} OBJECT ${CC_LIB_SRCS})
  else()
    add_library(${CC_LIB_NAME} INTERFACE ${CC_LIB_HDRS})
  endif()

  if(CC_LIB_DEPS)
    foreach(DEP IN LISTS CC_LIB_DEPS)
      target_link_libraries(${CC_LIB_NAME} PUBLIC ${DEP})
    endforeach()
  endif()
endfunction()

function(cc_test)
  set(options "")
  set(oneValueArgs NAME)
  set(multiValueArgs SRCS DEPS)
  cmake_parse_arguments(CC_LIB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN})
#  message("${CC_LIB_NAME}, ${CC_LIB_SRCS}, ${CC_LIB_DEPS}")
  set(EXEC "${CC_LIB_NAME}_EXEC")
  add_executable(${EXEC} ${CC_LIB_SRCS})

  target_link_libraries(${EXEC} PRIVATE ${GTEST_LIBRARIES} GTest::gtest_main Threads::Threads)
  foreach(DEP IN LISTS CC_LIB_DEPS)
    target_link_libraries(${EXEC} PRIVATE ${DEP})
  endforeach()

  add_test(NAME ${CC_LIB_NAME} COMMAND ${EXEC})
endfunction()