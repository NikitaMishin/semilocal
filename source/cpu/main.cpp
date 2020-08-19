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
//    std::cout<<sizeof(wordType) <<std::endl;
    std::srand(std::time(NULL)); // use current time as seed for random generator

//    std::vector<int> seq_a{1,0,0,1,1,0,1,1,0,0,1,0,0,1,1,0};
//    std::vector<int> seq_b{0,1,0,1,1,0,0,0 ,1,1,0,1,1,0,1,0};
//    ok a:00100010
//    ok b:01000101
//    ok a:1001101100100110
//    ok b:0101100011011010

// a:11111011100000000111010100011010
//b:10111000000000110000101101010011
//


    typedef size_t wordType;

    int thds = strtol(argv[1], NULL, 10);
//    int a_size = strtol(argv[2], NULL, 10);
//    int b_size = strtol(argv[3], NULL, 10);
    int a_size = 6400;
    int b_size = 6400*3;
//    10 10 10 00
//     11 01 00 00

//    64*640
//    64*6400
    auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({3,2,1,0}));
//    std::vector<int>seq_a = {1,2,3,4, 1,2,1,1};
//    std::vector<int>seq_b =  {2,2,2,2, 1,2,1,1};
//    std::vector<int>seq_b = {3,3,1,3,1,3,3,3};
//    std::vector<int>seq_a =  {1,1,1,1,1,1,1,1};

        auto seq_a = gen_vector_seq(a_size ,3);
    auto seq_b = gen_vector_seq( b_size,4);
    auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);

    auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
//
//    for (int i = 0; i < a_size; ++i) {
//        std::cout<<std::bitset<2>(mappers.first.at(seq_a[ seq_a.size() - 1 - i]))<<" ";
//    }
//    std::cout<<std::endl;
//    for (int i = 0; i < b_size; ++i) {
//        std::cout<<std::bitset<2> (mappers.first.at(seq_b[seq_b.size() - 1 - i]))<<" ";
//    }
//    std::cout<<std::endl;
//    std::cout<<std::endl;
////
//    std::cout<<std::bitset<8>(a.first.first[0])<<" "<<std::bitset<8>(a.first.first[1])<<std::endl;
//    std::cout<<std::bitset<8>(b.first.first[0])<<" "<<std::bitset<8>(b.first.first[1])<<std::endl;

int i = 0;
    while (true){
        i++;
        auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({3,2,1,0}));
        std::srand(NULL); // use current time as seed for random generator

        auto seq_a = gen_vector_seq(a_size ,4);
        auto seq_b = gen_vector_seq( b_size,4);

        auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);

        auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
        auto bit_mpi = prefix_lcs_via_braid_bits_4symbol_mpi(a.first.first,a.first.second, a.second ,
                                                    b.first.first,b.first.second, b.second,thds);
        auto c = prefix_lcs_sequential(seq_a,seq_b);

        if( c!=bit_mpi){

            std::cout<<"bit_mpi:"<<bit_mpi<<",bit:"<<",correct:"<<c<<std::endl;
            for (int j:seq_a) {
                std::cout<<j;
            }
            std::cout<<std::endl<<std::endl;
            for (int j:seq_b) {
                std::cout<<j;
            }
            std::cout<<std::endl;
            return 0 ;
        }

    }

    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << std::endl<<"Res4symbolalppabet: "<< prefix_lcs_via_braid_bits_4symbol_mpi(a.first.first,a.first.second, a.second ,
                                                      b.first.first,b.first.second, b.second,thds) << std::endl;
    auto time = std::chrono::high_resolution_clock::now() - begin;
    std::cout <<"Time 4symbol alphabet with open mpi:" <<std::chrono::duration<double, std::milli>(time).count() << std::endl;

//    std::cout << std::endl<<"Res binary mpi : "<< prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
////                                                                                   b.first.first,b.first.second, b.second,2,thds) << std::endl;


    mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({1,0}));
    a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);

    b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);


    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << std::endl<<": "<< prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
                                                                         b.first.first,b.first.second, b.second,thds) << std::endl;
    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
    std::cout <<"Time: binary:" <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;


auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << std::endl<<"res corret: "<< prefix_lcs_sequential(seq_a,seq_b) << std::endl;
    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
    std::cout <<"Time prefix_lcs:" <<std::chrono::duration<double, std::milli>(time2).count() << std::endl;






//
//
//    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"ResBin: "<< prefix_lcs_via_braid_bits_binary(a.first.first,a.first.second, a.second ,
//                                                      b.first.first,b.first.second, b.second,2) << std::endl;
//    auto time = std::chrono::high_resolution_clock::now() - begin;
//    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time).count() << std::endl;
//
//    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
////                                         b.first.first,b.first.second, b.second,2,thds)
//    std::cout <<"Res:"<< prefix_lcs_sequential(seq_a,seq_b) << std::endl;
//    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    std::cout << std::chrono::duration<double, std::milli>(time2).count() << std::endl;
//
//        auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << sticky_braid_mpi(seq_a, seq_b,1) << std::endl;
//    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//    std::cout << std::chrono::duration<double, std::milli>(time1).count() << std::endl;
////
//    auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << sticky_braid_sequential (seq_a, seq_b) << std::endl;
//    auto time3 = std::chrono::high_resolution_clock::now() - begin3;
//    std::cout << std::chrono::duration<double, std::milli>(time3).count() << std::endl;
//



    return 0 ;


}

double get_elapsed_time_ms() {
    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto time = std::chrono::high_resolution_clock::now() - begin;
    return std::chrono::duration<double, std::milli>(time).count();
}





