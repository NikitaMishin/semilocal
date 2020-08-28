cmake_minimum_required(VERSION 3.16)
project(cpu)

set(CMAKE_CXX_STANDARD 11)

add_library(cpu  library.h library_correctnesses_tests.h transposition_network_approach/transposition_network_4symbol_alphabet_bit.h fasta_parser.h library.cpp transposition_network_approach/encoders_and_decoders.h sequence_generators.h transposition_network_approach/transposition_network_binary_alphabet.h transposition_network_approach/transposition_network_unbounded_alphabet.h)
add_subdirectory(Google_tests)
#add_subdirectory(cpu_impl/clang10/prefix_lcs)

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
set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS}  -std=c++11 -fopenmp    -Wall -fstrict-aliasing   -ffast-math  -march=native  -O0    ")
#set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -std=c++11 -fopenmp -floop-parallelize-all -ftree-parallelize-loops=4    -Wall  -ffast-math  -march=native  -O3 -fstrict-aliasing  -fopt-info-vec-optimized-all -funsafe-math-optimizations     -fopt-info       ")
#-O3-g0 -Ofast  -ftree-vectorize -ffast-math -mprefer-vector-width=512 -march=skylake-avx512-g
#-mprefer-vector-width=512
#-fopenmp=libomp dd
#-vec -simd -qopenmp -O3
#-qopt-report=5
# -fopenmp-simd -simd- -fopt-info-vec-missed
#-fopt-info
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0")
add_executable(main main.cpp)
target_link_libraries(main cpu)


add_executable(prefix_lcs cpu_impl/clang10/prefix_lcs/prefix_lcs.cpp)
target_link_libraries(prefix_lcs cpu)


#CUDACXX=/usr/local/cuda-11.0/bin/nvcc
