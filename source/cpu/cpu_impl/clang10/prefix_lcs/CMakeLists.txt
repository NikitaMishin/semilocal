cmake_minimum_required(VERSION 3.16)

project(test_prefix_lcs)

set(CMAKE_CXX_STANDARD 11)



add_executable(test_prefix_lcs  prefix_lcs.cpp)
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -std=c++11 -fopenmp    -Wall -fstrict-aliasing   -ffast-math  -march=native  -O3    ")

target_link_libraries(test_prefix_lcs cpu )

