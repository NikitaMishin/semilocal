//
// Created by nikita on 27.07.2020.
//

#include "semi_local.h"
#include "transposition_network_approach/transposition_network_binary_alphabet.h"
#include "transposition_network_approach/encoders_and_decoders.h"
#include "sequence_generators.h"
#include "fasta_parser.h"
#include "naive_prefix_lcs.h"
#include "transposition_network_approach/4symbol_new.h"
#include "transposition_network_approach/2symbol_new.h"
#include "transposition_network_approach/2symbol_new_2.h"
#include <string>
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <algorithm>
#include <unordered_set>
//static const int length = 1024*8;
//static float a[length];





int main(int argc, char *argv[]) {

    typedef  unsigned int wordType;


//    std::cout<<sizeof(wordType) <<std::endl;
    std::srand(std::time(NULL)); // use current time as seed for random generator





    int a_size = 64*100*25;//
    int b_size = 64*100*26;
//    int a_size = 1024*32*10;//
//    int b_size = 1024*1024*64;
//96458

    auto seq_a = gen_vector_seq(a_size ,4);
    auto seq_b = gen_vector_seq(b_size,4);


    auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({3,2,1,0}));
    auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
    auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);

//    mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({1,0}));
//    a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
//    b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);


//


//
//    auto begin0 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"res: "<< prefix_lcs_via_braid_bits_2symbol_v3_full_mask<unsigned int,wordType>
//            (a.first.first,a.first.second, a.second ,b.first.first,b.first.second, b.second,1) << std::endl;
//    auto time0 = std::chrono::high_resolution_clock::now() - begin0;
//    std::cout <<"time 2: " <<std::chrono::duration<double, std::milli>(time0).count() << std::endl;


    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    std::cout << std::endl<<"res: "<<
              prefix_lcs_via_braid_bits_4symbol_v2_full_mask<wordType>(a.first.first,a.first.second, a.second ,
                                                                       b.first.first,b.first.second, b.second,1)<< std::endl;
    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
    std::cout <<"Time 4: " <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;

    std::cout << std::endl<<"res corret: "<< prefix_lcs_sequential(seq_a,seq_b) << std::endl;

    //
//
//    std::cout << std::endl<<"Res binary mpi : "<< prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
////                                                                                   b.first.first,b.first.second, b.second,2,thds) << std::endl;
//
//
//
//
//
//auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false

//    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    std::cout <<"Time prefix_lcs:" <<std::chrono::duration<double, std::milli>(time2).count() << std::endl;
//
//
//
//
//
//
////
////
//    auto begin4 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"ResBin: "<< prefix_lcs_via_braid_mpi_bitwise_operator<int,int>(seq_a,seq_b,1) << std::endl;
//    auto time4 = std::chrono::high_resolution_clock::now() - begin4;
//    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time4).count() << std::endl;
//
//    auto begin5 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"ResBin: "<< prefix_lcs_via_braid_mpi_less_operator <int,int>(seq_a,seq_b,1) << std::endl;
//    auto time5 = std::chrono::high_resolution_clock::now() - begin5;
//    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time5).count() << std::endl;
//
////
////    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//////    prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
//////                                         b.first.first,b.first.second, b.second,2,thds)
////    std::cout <<"Res:"<< prefix_lcs_sequential(seq_a,seq_b) << std::endl;
////    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
////    std::cout << std::chrono::duration<double, std::milli>(time2).count() << std::endl;
////
////        auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << sticky_braid_mpi(seq_a, seq_b,1) << std::endl;
////    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
////    std::cout << std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//////
////    auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << sticky_braid_sequential (seq_a, seq_b) << std::endl;
////    auto time3 = std::chrono::high_resolution_clock::now() - begin3;
////    std::cout << std::chrono::duration<double, std::milli>(time3).count() << std::endl;
////
//
//
//    //       ---
//    //     //   \\
//    //    //     \\
//    //   //=======\\
//    //  //         \\
//    //              ||
//    //
//
//
//    return 0 ;

}

double get_elapsed_time_ms() {
    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto time = std::chrono::high_resolution_clock::now() - begin;
    return std::chrono::duration<double, std::milli>(time).count();
}





