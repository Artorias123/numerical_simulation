cc_library(
  NAME Types
  HDRS Types.h
)

cc_library(
  NAME Compare
  HDRS Compare.h
)

cc_test(
  NAME UtilsTest
  SRCS Test/TypesTest.cpp
       Test/CompareTest.cpp
  DEPS Types
       Compare
)
