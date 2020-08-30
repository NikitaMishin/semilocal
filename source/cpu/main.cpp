//
// Created by nikita on 27.07.2020.
//

#include "semi_local.h"
#include "transposition_network_approach/transposition_network_binary_alphabet.h"
#include "transposition_network_approach/encoders_and_decoders.h"
#include "sequence_generators.h"
#include "transposition_network_approach/transposition_network_unbounded_alphabet.h"
#include "fasta_parser.h"
#include <string>
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <algorithm>
#include <unordered_set>
//static const int length = 1024*8;
//static float a[length];





int main(int argc, char *argv[]) {

    typedef unsigned int wordType;

    int thds = strtol(argv[1], NULL, 10);
    std::string a_filepath = std::string(argv[2]);
    std::string b_filepath = std::string(argv[3]);
    auto name_content_a = parse_input_file(a_filepath);
    auto name_content_b = parse_input_file(b_filepath);

    auto seq_a = transform_to_int_vector(name_content_a.second.second);
    auto seq_b = transform_to_int_vector(name_content_b.second.second);
    if(seq_a.size()>seq_b.size()) std::swap(seq_a,seq_b);

    auto beg_preprocess = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({3,2,1,0}));
    auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
    auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
    auto time_preprocess = std::chrono::high_resolution_clock::now() - beg_preprocess;
    auto elapsed_time_preprocess  = long (std::chrono::duration<double, std::milli>(time_preprocess).count());


    auto beg = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto res = prefix_lcs_via_braid_bits_4symbol(a.first.first,a.first.second, a.second ,b.first.first,b.first.second, b.second);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time  = long (std::chrono::duration<double, std::milli>(time).count());

    std::cout<<elapsed_time_preprocess<<std::endl;
    std::cout<<elapsed_time<<std::endl;
    std::cout<<res<<std::endl;
    std::cout<<seq_a.size()<<std::endl;
    std::cout<<seq_b.size();
    return 0;

//    std::cout<<sizeof(wordType) <<std::endl;
    std::srand(std::time(NULL)); // use current time as seed for random generator



//    ok a:00100010
//    ok b:01000101
//    ok a:1001101100100110
//    ok b:0101100011011010

// a:11111011100000000111010100011010
//b:10111000000000110000101101010011
//

//
//    typedef unsigned  long wordType;
//
////    int thds = strtol(argv[1], NULL, 10);
////    int a_size = strtol(argv[2], NULL, 10);
////    int b_size = strtol(argv[3], NULL, 10);
//    int a_size = 64*10000;//
//    int b_size = 64*10000;
////    int a_size = 64*100-31;
////    int b_size = 12573;
//
////    10 10 10 00
////     11 01 00 00
////    int a_size =4+1 ;//2320
////    int b_size = 4*3+4;
//
////    int a_size =64+1 ;//2320
////    int b_size = 32*4+32;
//
////    6539;//
////    int b_size = 6581;
//
////    64*640
////    64*6400
//    auto seq_a = gen_vector_seq(a_size ,2);
//    auto seq_b = gen_vector_seq( b_size,2);
//
//    auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({3,2,1,0}));
//    auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
//
//    auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
////
////    for (int i = 0; i < a_size; ++i) {
////        std::cout<<std::bitset<2>(mappers.first.at(seq_a[ seq_a.size() - 1 - i]))<<" ";
////    }
////    std::cout<<std::endl;
////    for (int i = 0; i < b_size; ++i) {
////        std::cout<<std::bitset<2> (mappers.first.at(seq_b[seq_b.size() - 1 - i]))<<" ";
////    }
////    std::cout<<std::endl;
////    std::cout<<std::endl;
//////
////    std::cout<<std::bitset<8>(a.first.first[0])<<" "<<std::bitset<8>(a.first.first[1])<<std::endl;
////    std::cout<<std::bitset<8>(b.first.first[0])<<" "<<std::bitset<8>(b.first.first[1])<<std::endl;
//
////int i = 0;
////    while (true){
////        std::cout<<"run number"<<i<<std::endl;
////        i++;
////        auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({3,2,1,0}));
////        std::srand(NULL); // use current time as seed for random generator
////
////        auto seq_a = gen_vector_seq(a_size ,2);
////        auto seq_b = gen_vector_seq( b_size,4);
////
////        auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
////
////        auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
////        auto bit_mpi = prefix_lcs_via_braid_bits_4symbol_mpi(a.first.first,a.first.second, a.second ,
////                                                    b.first.first,b.first.second, b.second,thds);
////        auto c = prefix_lcs_sequential(seq_a,seq_b);
////
////        if( c!=bit_mpi){
////
////            std::cout<<"bit_mpi:"<<bit_mpi<<",bit:"<<",correct:"<<c<<std::endl;
////            for (int j:seq_a) {
////                std::cout<<j;
////            }
////            std::cout<<std::endl<<std::endl;
////            for (int j:seq_b) {
////                std::cout<<j;
////            }
////            std::cout<<std::endl;
////            return 0 ;
////        }
////
////    }
//
//
//    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<": "<<
//                                prefix_lcs_via_braid_bits_4symbol(a.first.first,a.first.second, a.second ,
//                                                                 b.first.first,b.first.second, b.second)
//            << std::endl;
//    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//    std::cout <<"Time: 4symbol:" <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//
//    mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({1,0}));
//    a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
//
//    b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
//
//
//    auto begin = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"Binary: "<< prefix_lcs_via_braid_bits_binary(a.first.first,a.first.second, a.second ,
//                                                      b.first.first,b.first.second, b.second) << std::endl;
//    auto time = std::chrono::high_resolution_clock::now() - begin;
//    std::cout <<"Time binary:" <<std::chrono::duration<double, std::milli>(time).count() << std::endl;
//
//    auto begin0 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"Sequentual SImple: "<< prefix_lcs_via_braid_sequential<int,int>(seq_a,seq_b) << std::endl;
//
//    auto time0 = std::chrono::high_resolution_clock::now() - begin0;
//    std::cout <<"Simple:" <<std::chrono::duration<double, std::milli>(time0).count() << std::endl;
//
//
//
//
////    std::cout << std::endl<<"Res binary mpi : "<< prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
//////                                                                                   b.first.first,b.first.second, b.second,2,thds) << std::endl;
//
//
//
//
//
//auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"res corret: "<< prefix_lcs_sequential(seq_a,seq_b) << std::endl;
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





