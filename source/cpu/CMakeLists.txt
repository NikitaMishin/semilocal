cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_C_COMPILER /usr/bin/)
#set(CMAKE_CXX_COMPILER /usr/bin/g++-10)
#library_correctnesses_tests.h
add_library(cpu semi_local.h transposition_network_approach/transposition_network_4symbol_alphabet_bit.h fasta_parser.h semi_local.cpp transposition_network_approach/encoders_and_decoders.h sequence_generators.h transposition_network_approach/transposition_network_binary_alphabet.h transposition_network_approach/transposition_network_unbounded_alphabet.h transposition_network_approach/2symbol_new.h   transposition_network_approach/4symbol_new.h transposition_network_approach/2symbol_new_2.h unit_monge_mult/steady_ant.h)
add_subdirectory(Google_tests)

#/usr/bin/clang++-10
#/usr/bin/g++-10
#-ftree-vectorizer-verbose=
# added -fopenmp
# disable vectorization: -fno-tree-vectorize
# something to test: -fdisable-tree-cunrolli
# -fstrict-aliasing
#-fno-tree-vectorize
#-fno-tree-vectorize
#-fopenmp=libomp
#-fopt-info-vec-optimized-optimized
# -fopenmp=libiomp5
#-floop-parallelize-all
set(CMAKE_CXX_FLAGS   "${CMAKE_CXX_FLAGS}  -std=c++11 -fopenmp  -Wall     -fstrict-aliasing   -fopt-info-vec-optimized-optimized -ffast-math -fopenmp   -march=native  -O0  -fopt-info-vec-optimized-all -funsafe-math-optimizations -fopt-info-vec-optimized-optimized " )#    ")
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -floop-parallelize-all -ftree-parallelize-loops=4    -Wall  -ffast-math  -march=native  -O3 -fstrict-aliasing  -fopt-info-vec-optimized-all -funsafe-math-optimizations     -fopt-info       ")
#-O3-g0 -Ofast  -ftree-vectorize -ffast-math -mprefer-vector-width=512 -march=skylake-avx512-g
#-mprefer-vector-width=512
#-fopenmp=libomp dd
#-vec -simd -qopenmp -O3
#-qopt-report=5
# -fopenmp-simd -simd- -fopt-info-vec-missed
#-fopt-info
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_C XX_FLAGS "${CMAKE_CXX_FLAGS} -O0")


add_executable(main main.cpp)
target_link_libraries(main cpu)


#library.h transposition_network_approach/transposition_network_4symbol_alphabet_bit.h fasta_parser.h library.cpp transposition_network_approach/encoders_and_decoders.h sequence_generators.h transposition_network_approach/transposition_network_binary_alphabet.h transposition_network_approach/transposition_network_unbounded_alphabet.h
# prefix lcs sequential
add_executable(prefix_lcs cpu_impl/prefix_lcs.cpp )

add_executable(semi_local_lcs_sequential cpu_impl/semi_local_lcs_sequential.cpp )

add_executable(semi_local_lcs_mpi cpu_impl/semi_local_lcs_mpi.cpp )

add_executable(prefix_lcs_via_braid_mpi_less_operator cpu_impl/prefix_lcs_via_braid_mpi_less_operator.cpp)

add_executable(prefix_lcs_via_braid_bits_4symbol_int cpu_impl/prefix_lcs_via_braid_bits_4symbol_int.cpp)
#CUDACXX=/usr/local/cuda-11.0/bin/nvcc

add_executable(prefix_lcs_via_braid_bits_4symbol_long_long cpu_impl/prefix_lcs_via_braid_bits_4symbol_long_long.cpp)

add_executable(prefix_lcs_via_braid_bits_4symbol_mpi_int cpu_impl/prefix_lcs_via_braid_bits_4symbol_mpi_int.cpp)

add_executable(prefix_lcs_via_braid_bits_4symbol_mpi_long_long cpu_impl/prefix_lcs_via_braid_bits_4symbol_mpi_long_long.cpp)

add_executable(prefix_lcs_via_braid_bits_binary_int cpu_impl/prefix_lcs_via_braid_bits_binary_int.cpp)

add_executable(prefix_lcs_via_braid_bits_binary_long_long cpu_impl/prefix_lcs_via_braid_bits_binary_long_long.cpp)

add_executable(prefix_lcs_via_braid_bits_binary_mpi_int cpu_impl/prefix_lcs_via_braid_bits_binary_mpi_int.cpp)

add_executable(prefix_lcs_via_braid_bits_binary_mpi_long_long cpu_impl/prefix_lcs_via_braid_bits_binary_mpi_long_long.cpp)

add_executable(prefix_lcs_via_braid_sequential cpu_impl/prefix_lcs_via_braid_sequential.cpp)