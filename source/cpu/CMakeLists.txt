cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)

add_library(cpu  library.h  prefix_lcs.h transposition_network.h utils.h library.cpp)

# added -fopenmp
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp   -Wall -Ofast  ")
#-O3-g0

#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
add_executable(main main.cpp)
target_link_libraries(main cpu)

#CUDACXX=/usr/local/cuda-11.0/bin/nvcc