project(UnscentedKF)

cmake_minimum_required (VERSION 3.5)

add_definitions(-std=c++11)

add_subdirectory(
        lib/googletest)

include_directories(
        lib/googletest/include)

set(lib_sources
        src/ukf.cpp
        src/ukf.h
        src/tools.cpp
        src/tools.h)

set(lib_test_sources
        src/ukf_test.cpp
        src/tools_test.cpp)

add_library(UnscentedKFLib ${lib_sources})
add_executable(UnscentedKF src/main.cpp)
add_executable(UnscentedKFTests ${lib_test_sources})
target_link_libraries(UnscentedKF UnscentedKFLib)
target_link_libraries(UnscentedKFTests UnscentedKFLib gtest gtest_main)