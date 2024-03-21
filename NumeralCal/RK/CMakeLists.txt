cc_library(
  NAME RK
  HDRS RK.h
)

cc_test(
  NAME RKTest
  SRCS Test/RKTest.cpp
  DEPS RK
)

# add_executable(RKTest Test/RKTest.cpp)
# target_link_libraries(RKTest ${GTEST_LIBRARIES} gtest_main pthread)
# 
# add_test(NAME test_rk COMMAND RKTest)

# add_library(Engine SHARED ${EngineSrc})
# target_compile_options(Engine PRIVATE -fvisibility=hidden)
# target_link_libraries(Engine PRIVATE Item)
