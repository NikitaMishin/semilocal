

#include <chrono>
#include <cstdio>
#include <iostream>
#include "memory_management.h"
#include "cpu/semi_local.h"

//#include "first_iteration/transposition_network_4symbol_gpu.cu"
//#include "first_iteration/semi_local.cu"


//#include "kawanami.cuh"
//#include "types.h"
//#include "utils/transformers.h"
//#include "semi_local/algorithms.cuh"
//using namespace memory_management;

#include "cpu_routines/cpu_routines.h"
#include "gpu_routines/gpu_routines.cuh"

#include "semi_local/algorithms.cuh"
#include "utils/sequence_generators.h"


int main() {
    int a_size = 1305;
    int b_size = 1000;

    auto seq_a = gen_array(a_size,23);
    auto seq_b = gen_array(b_size,23);


    auto strategy = SemiLocalAntiDiagonalStrategy<int,true,false,false>(4,128,400,300);
    auto should = Permutation(a_size + b_size, a_size + b_size);
    auto actual = Permutation(a_size + b_size, a_size + b_size);

    strategy.compute(actual,seq_a,a_size,seq_b,b_size);

    semi_local::sticky_braid_mpi<int,true,true>(should,seq_a,a_size,seq_b,b_size,4);

    std::cout<<should.is_equal_to(actual);


//    typedef unsigned long long wordType ;
//
////
////    int a_size = 64*100;//
////    int b_size = 64*100;
//
//    int a_size = 10240000;//
//    int b_size = 10240000;
//
//    a_size = a_size;//
//    b_size = b_size;
//
//    //bitsetleft
//
//    auto bitset_left = allocate_1D_array_aligned_on_gpu<wordType>(a_size+b_size);
////    auto bitset_top = allocate_1D_array_aligned_on_gpu<wordType>(b_size);
//    auto a_rev = allocate_1D_array_aligned_on_gpu<wordType>(a_size);
//    auto b = allocate_1D_array_aligned_on_gpu<wordType>(b_size);
//
//
//    auto c = allocate_1D_array_aligned_on_cpu<wordType>(a_size);
//    auto d = allocate_1D_array_aligned_on_cpu<wordType>(b_size);
//    for(int i=0; i <a_size;i++){
//        c[i] = rand() % 40000;
//    }
//    for(int i=0; i <b_size;i++){
//        d[i] = rand() % 40000000;
//    }
//    copy_from_cpu_to_gpu_sync(c,a_rev,a_size);
//    copy_from_cpu_to_gpu_sync(d,b,b_size);
//
//    for(int i=0; i <a_size;i++){
//        c[i] = rand() % 4000000;
//    }
//    for(int i=0; i <b_size;i++){
//        d[i] = rand() % 40000000;
//    }
//    copy_from_cpu_to_gpu_sync(c,bitset_left,a_size);
////    copy_from_cpu_to_gpu_sync(d,bitset_top,b_size);
//
//    auto begin0 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    semi_local_lcs_gpu::gpu_wrappers::semi_local_process_diagonal(a_rev,a_size,b,b_size,bitset_left,512,256*100,0,0,a_size);
//    auto time0 = std::chrono::high_resolution_clock::now() - begin0;
//    memory_management::gpuAssert(cudaGetLastError(), __FILE__, __LINE__);
////    copy_from_gpu_to_cpu_sync(res,c,a_size);
//    copy_from_gpu_to_cpu_sync(bitset_left,c,a_size);
//
//    std::cout<<c[0]<<std::endl;
//    std::cout <<"time 2: " <<std::chrono::duration<double, std::milli>(time0).count() << std::endl;
//

    return 0;
}