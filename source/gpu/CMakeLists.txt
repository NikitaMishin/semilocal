cmake_minimum_required(VERSION 3.16)
project(gpu CUDA)

#set(CMAKE_CUDA_COMPILER /usr/local/cuda-11.1/bin/nvcc)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --ptxas-options=-v -Xptxas -O3 -v")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_CXX_STANDARD 11)

string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_50,code=sm_50")


add_executable( gpu main.cu )
#CUDACXX=/usr/local/cuda-11.0/bin/nvcc

#/usr/local/cuda-11.0/bin/nvcc
set_target_properties(
        gpu
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)