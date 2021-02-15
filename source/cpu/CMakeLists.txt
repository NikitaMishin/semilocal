cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)

add_library(cpu semi_local.h  semi_local.cpp prefix_lcs/bitwise/transposition_network_binary_alphabet.h )
#add_subdirectory(Google_tests)

# disable vectorization: -fno-tree-vectorize

#set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++11 -O0 -fopenmp  -fno-tree-vectorize -fopt-info-vec-optimized-optimized -fopt-info-vec-optimized-all")
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++11 -fopenmp  -Wall       -march=native  -O3  " )#    ")

add_executable(main main.cpp)
target_link_libraries(main cpu)

##
add_executable(prefix_lcs cpu_impl/prefix_lcs.cpp)

##
add_executable(prefix_lcs_skewed cpu_impl/prefix_lcs_skewed.cpp)
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
