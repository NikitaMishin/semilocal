cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)

add_library(cpu semi_local.h
        semi_local.cpp prefix_lcs/bitwise/transposition_network_binary_alphabet.h
        prefix_lcs/bitwise/transposition_network_4symbol_alphabet_bit.h
        unit_monge_mult/smawk.h approximate_matching/CompleteAMatch.h)
#add_subdirectory(Google_tests)

# disable vectorization: -fno-tree-vectorize

#set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++11 -O0 -fopenmp  -fno-tree-vectorize -fopt-info-vec-optimized-optimized -fopt-info-vec-optimized-all")
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++17 -fopenmp  -march=native  -O3  " )#    ")
# -I /home/garrancha/libs/boost_1_75_0/ -for  boost libary
#add_executable(main main.cpp)
#target_link_libraries(main cpu)

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
add_executable(semi_local_parallel_withif_iterative cpu_impl/semi_local_parallel_withif_iterative.cpp)


#
##
add_executable(semi_local_parallel_1and3_combined_iterative cpu_impl/semi_local_parallel_1and3_combined_iterative.cpp)
#
##
#add_executable(semi_local_naive_rec cpu_impl/semi_local_naive_rec.cpp)
#
#
#
##
add_executable(semi_local_parallel_hybrid cpu_impl/semi_local_parallel_hybrid.cpp)
#

add_executable(semi_local_hybrid_iterative cpu_impl/semi_local_hybrid_iterative.cpp)
#




#
##

add_executable(prefix_lcs_semi_bit_parallel_1form cpu_impl/prefix_lcs_semi_bit_parallel_1form.cpp)


##
add_executable(prefix_lcs_semi_bit_parallel_2form cpu_impl/prefix_lcs_semi_bit_parallel_2form.cpp)

##
add_executable(prefix_lcs_semi_bit_parallel_old cpu_impl/prefix_lcs_semi_bit_parallel_old.cpp)


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

find_package(GTest REQUIRED )
include_directories(${GTEST_INCLUDE_DIRS})

add_executable(
        ApproximateMatchingTest
        test/ApproximateMatchingTest.cpp
)
target_link_libraries(ApproximateMatchingTest GTest::gtest GTest::gtest_main  cpu  )
target_include_directories(ApproximateMatchingTest  PUBLIC cpu ${GTEST_INCLUDE_DIRS}  )

