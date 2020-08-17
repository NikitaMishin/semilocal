//
// Created by nikita on 27.07.2020.
//

#include "library.h"
#include <string>
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <algorithm>
//static const int length = 1024*8;
//static float a[length];




int main(int argc, char *argv[]) {

    std::srand(std::time(NULL)); // use current time as seed for random generator
    int thds = strtol(argv[1], NULL, 10);
    int a_size = strtol(argv[2], NULL, 10);
    int b_size = strtol(argv[3], NULL, 10);
//    64*640
//    64*6400
    auto mappers = encode_alphabet<int,size_t>(std::unordered_set<int>({1,0}));
    auto seq_a = gen_vector_seq(a_size ,2);
    auto seq_b = gen_vector_seq( b_size,2);
//    std::vector<int> seq_a{1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0};
//    std::vector<int> seq_b{0,1,0,1,1,0,0,0,1,1,0,1,1,0,1,0};
//    ok a:00100010
//    ok b:01000101
//    ok a:1001101100100110
//    ok b:0101100011011010

// a:11111011100000000111010100011010
//b:10111000000000110000101101010011
//

std::cout<<thds<<std::endl;
    auto a = encode_reverse<int,size_t>(seq_a,&mappers.first,&mappers.second);

    auto b = encode<int, size_t >(seq_b,&mappers.first,&mappers.second);



//    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"ResBin: "<< prefix_lcs_via_braid_bits_binary(a.first.first,a.first.second, a.second ,
//                                                      b.first.first,b.first.second, b.second,2) << std::endl;
//    auto time = std::chrono::high_resolution_clock::now() - begin;
//    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time).count() << std::endl;

    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout <<"Res:"<< prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
                                                          b.first.first,b.first.second, b.second,2,thds) << std::endl;
    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    auto end = omp_get_wtime();
//    printf("Time: %f\n", end - start);
//    std::cout<<std::endl;
    std::cout << std::chrono::duration<double, std::milli>(time2).count() << std::endl;

        auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << prefix_lcs_sequential(seq_a, seq_b) << std::endl;
    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
    std::cout << std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//
//    auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << prefix_lcs_sequential(seq_a, seq_b) << std::endl;
//    auto time3 = std::chrono::high_resolution_clock::now() - begin3;
//    std::cout << std::chrono::duration<double, std::milli>(time3).count() << std::endl;




    return 0 ;


}

double get_elapsed_time_ms() {
    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto time = std::chrono::high_resolution_clock::now() - begin;
    return std::chrono::duration<double, std::milli>(time).count();
}





