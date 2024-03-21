cc_library(
  NAME RK
  HDRS RK.h
       ButcherTable.h
)

cc_test(
  NAME RKTest
  SRCS Test/RKTest.cpp
       Test/ButcherTableTest.cpp
  DEPS RK
)
