//
// Created by nikita on 27.07.2020.
//

#include "library.h"
#include <string>
#include <iostream>
#include <chrono>
#include <immintrin.h>


#include "omp.h"

inline __m256d reverse(__m256d x) {
    x = _mm256_permute2f128_pd(x, x, 1);
    x = _mm256_permute_pd(x, 5);
    return x;
}

int main() {

//    #pragma omp parallel
//    {
//        int ID = omp_get_thread_num();
//        std::cout<<"Hello from id="<<ID<<std::endl;
//        std::cout<<omp_get_num_procs();
//
//    }
//    std::srand(std::time(nullptr)); // use current time as seed for random generator
    std::srand(std::time(0)); // use current time as seed for random generator

    std::string str_a = "10", str_b = "010";

    auto a = gen_vector_seq(5000 * 2, 4345435);
    auto b = gen_vector_seq(100000 * 5, 3435435345);


    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << sticky_braid_sequential(a, b) << std::endl;
    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
    std::cout << std::chrono::duration<double, std::milli>(time2).count() << std::endl;


    auto cor = sticky_braid_sequential(a, b);
    auto cur = sticky_braid_sequential_skewed(a, b);

//    for (int i = 0; i <a.size()+b.size() ; ++i) {
//
//        if (cor[i]!=cur[i]) {
//            std::cout<<"Bad in "<< i <<std::endl;
//            return 1;
//        }
//
//    }

//
//
//
    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << sticky_braid_sequential_skewed(a, b) << std::endl;
    auto time = std::chrono::high_resolution_clock::now() - begin;
    std::cout << std::chrono::duration<double, std::milli>(time).count() << std::endl;

    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << prefix_lcs_sequential(a, b) << std::endl;
    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
    std::cout << std::chrono::duration<double, std::milli>(time1).count() << std::endl;



//    __m256d x = _mm256_set_pd(13,12,11,10);
//
//    std::cout << x.m256d_f64[0] << "  " << x.m256d_f64[1] << "  " << x.m256d_f64[2] << "  " << x.m256d_f64[3] <<std::endl;
//    x = reverse(x);
//    std::cout << x.m256d_f64[0] << "  " << x.m256d_f64[1] << "  " << x.m256d_f64[2] << "  " << x.m256d_f64[3] << std::endl;
//}

    return 0;
}

double get_elapsed_time_ms() {
    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto time = std::chrono::high_resolution_clock::now() - begin;
    return std::chrono::duration<double, std::milli>(time).count();
}





