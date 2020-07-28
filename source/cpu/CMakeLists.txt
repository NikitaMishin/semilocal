cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)

add_library(cpu library.cpp library.h)

add_executable(main main.cpp)
target_link_libraries(main cpu)