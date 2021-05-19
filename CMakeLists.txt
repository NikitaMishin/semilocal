cmake_minimum_required(VERSION 3.10.2)
project(semi CUDA)

set(APP_ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# set your 3 other root dirs, i.e. foo, main and test under app/src.
set(APP_ROOT_SRC_DIR ${APP_ROOT_DIR}/source)
set(APP_ROOT_COMMON_DIR ${APP_ROOT_SRC_DIR}/common)
set(APP_ROOT_COMMON_MATRICES_DIR ${APP_ROOT_COMMON_DIR}/monge)


set(SHARED_INCLUDES

#        ${APP_ROOT_COMMON_MATRICES_DIR}




        ${APP_ROOT_COMMON_DIR}
        #
#        ${APP_ROOT_COMMON_DIR}/monge




        # ${APP_ROOT_TEST_DIR}/<your-other-child-dirs>
        )

# This function will have effect to all the downstream cmakelist files.
include_directories(${SHARED_INCLUDES})




#set(CMAKE_CUDA_COMPILER /usr/local/cuda-11.1/bin/nvcc)
set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --ptxas-options=-v -Xptxas -O3 -v")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

set(CMAKE_CXX_STANDARD 11)

string(APPEND CMAKE_CUDA_FLAGS " -gencode arch=compute_50,code=sm_50")
#
#include_directories(${CMAKE_SOURCE_DIR}/source/gpu)


add_executable(semi source/gpu/runner.cu )
#CUDACXX=/usr/local/cuda-11.0/bin/nvcc

#/usr/local/cuda-11.0/bin/nvcc
set_target_properties(
        semi
        PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)