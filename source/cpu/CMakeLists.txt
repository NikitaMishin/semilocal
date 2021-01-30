cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_C_COMPILER /usr/bin/)
#set(CMAKE_CXX_COMPILER /usr/bin/g++-10)

add_library(cpu semi_local.h  semi_local.cpp )
#add_subdirectory(Google_tests)

# disable vectorization: -fno-tree-vectorize

#set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++11 -O0 -fopenmp  -fno-tree-vectorize -fopt-info-vec-optimized-optimized -fopt-info-vec-optimized-all")
#set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++11 -O1")
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++11 -fopenmp  -Wall    -fstrict-aliasing   -fopt-info-vec-optimized-optimized -ffast-math -fopenmp   -march=native  -O3  -fopt-info-vec-optimized-all -funsafe-math-optimizations -fopt-info-vec-optimized-optimized " )#    ")
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11  -march=skylake ")

add_executable(main main.cpp)
target_link_libraries(main cpu)

##
add_executable(prefix_lcs cpu_impl/prefix_lcs.cpp)
#
##
add_executable(semi_local_naive_iterative cpu_impl/semi_local_naive_iterative.cpp)
#
##
add_executable(semi_local_parallel_iterative cpu_impl/semi_local_parallel_iterative.cpp)
#
add_executable(semi_local_parallel_withoutif_iterative cpu_impl/semi_local_parallel_withoutif_iterative.cpp)


#
##
add_executable(semi_local_parallel_1and3_combined_iterative cpu_impl/semi_local_parallel_1and3_combined_iterative.cpp)
#
##
add_executable(semi_local_naive_rec cpu_impl/semi_local_naive_rec.cpp)
#
#
#
##
add_executable(semi_local_parallel_hybrid cpu_impl/semi_local_parallel_hybrid.cpp)
#


##
add_executable(braid_multiplication_parallel cpu_impl/braid_multiplication_parallel.cpp)
#
##
add_executable(braid_multiplication_sequential_non_optimized cpu_impl/braid_multiplication_sequential_non_optimized.cpp)
#
#
##
add_executable(braid_multiplication_sequential_memory cpu_impl/braid_multiplication_sequential_memory.cpp)
#
##
add_executable(braid_multiplication_sequential_precompute cpu_impl/braid_multiplication_sequential_precompute.cpp)
#
