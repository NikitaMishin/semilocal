cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)

add_library(cpu library.cpp library.h)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -O3  -Wall ")
#-Ofast -g0
add_executable(main main.cpp)
target_link_libraries(main cpu)

#CUDACXX=/usr/local/cuda-11.0/bin/nvcc