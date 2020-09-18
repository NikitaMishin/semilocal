cmake_minimum_required(VERSION 3.16)
project(gpu CUDA)

set(CMAKE_CUDA_COMPILER /usr/local/cuda-11.0/bin/nvcc)
set(CMAKE_CUDA_STANDARD 11)

add_executable( gpu main.cu first_iteration/semi_local.cu first_iteration/transposition_network_4symbol_gpu.cu)
#CUDACXX=/usr/local/cuda-11.0/bin/nvcc

#/usr/local/cuda-11.0/bin/nvcc
set_target_properties(
        gpu
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)