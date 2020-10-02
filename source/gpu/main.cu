

#include <chrono>
#include "memory_management.h"
#include "first_iteration/transposition_network_4symbol_gpu.cu"
using namespace memory_management;

int main() {

    typedef unsigned long long wordType ;

//
//    int a_size = 64*100;//
//    int b_size = 64*100;

    int a_size = 1024*32*10;//
    int b_size = 1024*1024*64;

    a_size = a_size / 32;//
    b_size = b_size / 32;

    //bitsetleft

    auto bitset_left = allocate_1D_array_aligned_on_gpu<wordType>(a_size);
    auto bitset_top = allocate_1D_array_aligned_on_gpu<wordType>(b_size);
    auto a_rev = allocate_1D_array_aligned_on_gpu<wordType>(a_size);
    auto b = allocate_1D_array_aligned_on_gpu<wordType>(b_size);


    auto c = allocate_1D_array_aligned_on_cpu<wordType>(a_size);
    auto d = allocate_1D_array_aligned_on_cpu<wordType>(b_size);
    for(int i=0; i <a_size;i++){
        c[i] = rand() % 40000;
    }
    for(int i=0; i <b_size;i++){
        d[i] = rand() % 40000000;
    }
    copy_from_cpu_to_gpu_sync(c,a_rev,a_size);
    copy_from_cpu_to_gpu_sync(d,b,b_size);

    for(int i=0; i <a_size;i++){
        c[i] = rand() % 4000000;
    }
    for(int i=0; i <b_size;i++){
        d[i] = rand() % 40000000;
    }
    copy_from_cpu_to_gpu_sync(c,bitset_left,a_size);
    copy_from_cpu_to_gpu_sync(d,bitset_top,b_size);

    auto begin0 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    four_symbol_gpu_runner_fully_gpu(a_rev,a_size,a_size*32,b,b_size,b_size*32,bitset_left,bitset_top,512);
    auto time0 = std::chrono::high_resolution_clock::now() - begin0;
    memory_management::gpuAssert(cudaGetLastError(), __FILE__, __LINE__);
//    copy_from_gpu_to_cpu_sync(res,c,a_size);
    copy_from_gpu_to_cpu_sync(bitset_left,c,a_size);

    std::cout<<c[0]<<std::endl;
    std::cout <<"time 2: " <<std::chrono::duration<double, std::milli>(time0).count() << std::endl;


    return 0;
}