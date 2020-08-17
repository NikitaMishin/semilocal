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
//    std::cout<<~0U<<std::endl;
////    std::printf("8-byte aligned address: %p\n", static_cast<void*>(encode<int>(gen_vector_ACGT(17)).first));
//    auto  seq = gen_vector_seq(15,2);
//
////      std::cout<<f_call.first.second.first<<" "<<f_call.first.second.second<<std::endl;
//
//    std::cout<<std::endl;
//    for (int i: seq) {
//        std::cout<<i;
//    }
//    std::cout<<std::endl;
//    auto s = encode_reverse<int,size_t>(seq,&mappers.first,&mappers.second);
//    auto d = decode_reverse(s.first.first,s.first.second,s.second,mappers.second);
//
//    for (int i: d) {
//        std::cout<<i;
//    }
//
////        auto a1 = gen_vector_seq( 100000/2, 40000000);
////    auto b1 = gen_vector_seq(500000*2/2, 200000000);
//////
////auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << prefix_lcs_sequential(a1, b1) << std::endl;
////    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
////    std::cout <<"Hahah"<< std::chrono::duration<double, std::milli>(time1).count() << std::endl;
////
    std::srand(std::time(NULL)); // use current time as seed for random generator

    auto mappers = encode_alphabet<int,unsigned long long>(std::unordered_set<int>({1,0}));
    auto seq_a = gen_vector_seq( 64*19,2);
    auto seq_b = gen_vector_seq( 64*20,2);
//    std::vector<int> seq_a{1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0};
//    std::vector<int> seq_b{0,1,0,1,1,0,0,0,1,1,0,1,1,0,1,0};
//    ok a:00100010
//    ok b:01000101
//    ok a:1001101100100110
//    ok b:0101100011011010

// a:11111011100000000111010100011010
//b:10111000000000110000101101010011
//

    auto a = encode_reverse<int,unsigned long long>(seq_a,&mappers.first,&mappers.second);

    auto b = encode<int, unsigned long long >(seq_b,&mappers.first,&mappers.second);

//    std::reverse(seq_a.begin(),seq_a.end());
//    std::reverse(seq_b.begin(),seq_b.end());
    std::cout<<"a:";
    for(auto elem : seq_a){
        std::cout<<elem;
    }
    std::cout<<std::endl<<"b:";
    for(auto elem : seq_b){
        std::cout<<elem;
    }


    std::cout<<std::endl;
    std::cout<<"a_rev:"<<std::bitset<sizeof(char)*8>(a.first.first[0]);
    std::cout<<std::bitset<sizeof(char)*8>(a.first.first[1])<<std::endl;
    std::cout<<"b:"<<std::bitset<sizeof(char)*8>(b.first.first[0]);
    std::cout<<std::bitset<sizeof(char)*8>(b.first.first[1])<<std::endl;

    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << std::endl<<"ResBin: "<< prefix_lcs_via_braid_bits_binary(a.first.first,a.first.second, a.second ,
                                                      b.first.first,b.first.second, b.second,2) << std::endl;
    auto time = std::chrono::high_resolution_clock::now() - begin;
    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time).count() << std::endl;

    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout <<"Res:"<< prefix_lcs_sequential(seq_a, seq_b) << std::endl;
    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    auto end = omp_get_wtime();
//    printf("Time: %f\n", end - start);
//    std::cout<<std::endl;
    std::cout << std::chrono::duration<double, std::milli>(time2).count() << std::endl;

    return 0 ;
//    std::cout<<sizeof(size_t)<<std::endl;
    int threads = strtol(argv[1], NULL, 10);
//    std::cout<<threads;
//    #pragma omp parallel
//    {
//        int ID = omp_get_thread_num();
//        std::cout<<"Hello from id="<<ID<<std::endl;
//        std::cout<<omp_get_num_procs();
//
//    }
//    std::srand(std::time(nullptr)); // use current time as seed for random generator

//    std::string str_a = "10", str_b = "010";
//    auto a = gen_vector_seq( 100000/2, 40000000);
//    auto b = gen_vector_seq(500000*2, 200000000);


    //    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << prefix_lcs_sequential(a, b) << std::endl;
//    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//    std::cout << std::chrono::duration<double, std::milli>(time1).count() << std::endl;


//    auto start = omp_get_wtime();
//    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << sticky_braid_mpi(a, b,threads) << std::endl;
//    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    auto end = omp_get_wtime();
//    printf("Time: %f\n", end - start);
//    std::cout<<std::endl;
//    std::cout << std::chrono::duration<double, std::milli>(time2).count() << std::endl;

//std::cout<<std::endl;
//    auto cor = sticky_braid_sequential(a, b);
//    auto cur = sticky_braid_mpi(a, b,8);
//    for (int i = 0; i <a.size()+b.size(); ++i) {
//        if (cor[i]!=cur[i])std::cout<<"ddddd";
//    }


//
//
//
//    std::reverse(a.begin(),b.end());



    return 0;
}

double get_elapsed_time_ms() {
    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto time = std::chrono::high_resolution_clock::now() - begin;
    return std::chrono::duration<double, std::milli>(time).count();
}





