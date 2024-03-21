cmake_minimum_required(VERSION 3.15)
project(game)

# 设置C++标准
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
add_compile_options(-fPIC)

# 设置生成 compile_commands.json
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})
enable_testing()

# 配置clang-tidy
# set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
add_custom_target(
    clang-tidy
    COMMAND run-clang-tidy -fix -format -quiet -j 4
    COMMENT "Running clang-tidy"
)

# 配置iwyu
find_program(IWYU_PATH NAMES include-what-you-use include-what-you-use)
add_custom_target(
    iwyu
    COMMAND python3 ${CMAKE_SOURCE_DIR}/build_tools/iwyu_tool.py -p . -- -Xiwyu --transitive_includes_only -Xiwyu --cxx17ns -Xiwyu --mapping_file=${CMAKE_SOURCE_DIR}/build_tools/iwyu.imp
    COMMENT "Running iwyu"
)

add_custom_target(
    clang-format
    COMMAND bash ${CMAKE_SOURCE_DIR}/build_tools/clang_format.sh ${CMAKE_SOURCE_DIR}
    COMMENT "Running clang-format"
)

# 添加工程的根目录到头文件搜索路径
include_directories(${CMAKE_SOURCE_DIR})
include(${CMAKE_CURRENT_SOURCE_DIR}/build_tools/macros.cmake)

add_subdirectory(NumeralCal)
