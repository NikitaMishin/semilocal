//
// Created by nikita on 27.07.2020.
//
//
#include "semi_local.h"
#include "prefix_lcs/bitwise/encoders_and_decoders.h"
#include "sequence_generators.h"
#include "naive_prefix_lcs.h"
#include "prefix_lcs/bitwise/4symbol_new.h"
#include "prefix_lcs/bitwise/2symbol_new.h"
#include "prefix_lcs/bitwise/2symbol_new_2.h"
#include "unit_monge_mult/steady_ant.h"
#include "prefix_lcs/bitwise/transposition_network_binary_alphabet.h"
#include "prefix_lcs/bitwise/transposition_network_4symbol_alphabet_bit.h"
#include <string>
#include <iostream>
#include <chrono>
#include <immintrin.h>
#include <algorithm>
#include <unordered_set>
//static const int length = 1024*8;
//static float a[length];
//
//
//#include "boost/multiprecision/cpp_int.hpp"

int f() { return 5; };

int main() {

    typedef unsigned short Hold;
    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
    precalc(map, 5);

    omp_set_nested(true);
    auto a_size = 70;
    auto b_size = 125;

    auto seq_a = gen_vector_seq<Hold>(a_size, 26);
    auto seq_b = gen_vector_seq<Hold>(b_size, 26);
    auto a = new Hold[a_size];
    auto b = new Hold[b_size];
    auto a_int = new int[a_size];
    auto b_int = new int[b_size];
    for (int i = 0; i < a_size; ++i) {
        a[i] = seq_a[i];
        a_int[i] = seq_a[i];
//        std::cout<<a[i]<<" ";
    }
    for (int i = 0; i < b_size; ++i) {
        b[i] = seq_b[i];
        b_int[i] = seq_b[i];
//        std::cout<<b[i]<<" ";
    }

//6555
//    time24866.5
//    Precalc for value 5
//    time24010.9
//    1
// 26
// 40
    auto should = Permutation(a_size + b_size, a_size + b_size);
    auto actual = Permutation(a_size + b_size, a_size + b_size);

//

//    std::cout << "Precalc for value 5" << std::endl;
//    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    semi_local::hybrid_iterative<Hold , false> (actual, a, a_size, b, b_size, map,5,20,4);
//    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//    std::cout << "total_time" << std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//
//
//        std::cout << "Precalc for value 5" << std::endl;
//    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    semi_local::sticky_braid_mpi<int , false, true>(should, a_int, a_size, b_int, b_size, 4);
//
//    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    std::cout << "time" << std::chrono::duration<double, std::milli>(time2).count() << std::endl;


//    std::cout << should.is_equal_to(actual);
//

    std::cout<<prefix_lcs_sequential(a_int,a_size,b_int,b_size)<<std::endl;
    std::cout<<prefix_lcs_sequential_skewed(a_int,a_size,b_int,b_size)<<std::endl;
    delete []a;
    delete [] b;
    delete [] a_int;
    delete [] b_int;
//    int a[a_size];
//    int b[b_size];
//    for (int i = 0; i < a_size; ++i) {
//        a[i] = seq_a[i];
////        std::cout<<a[i]<<" ";
//    }
//    for (int i = 0; i < b_size; ++i) {
//        b[i] = seq_b[i];
////        std::cout<<b[i]<<" ";
//    }
//    std::cout<<std::endl;
//
//    //
//    // |\ /|
//
//
//

//
//
//
//    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//
//    auto should = Permutation(a_size+b_size,a_size+b_size);
//
//
//
//
//    semi_local::hybrid_approach::first_and_third_phase_merged(a, a_size, b, b_size, should, map, 2,32);
//    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    std::cout << "merged " << std::chrono::duration<double, std::milli>(time2).count()
//              << std::endl;
//
//
//    auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//
//    auto actual = new  Permutation(a_size+b_size,a_size+b_size);
//    semi_local::strand_combing_approach::sticky_braid_mpi(*actual,a,a_size,b,b_size,32);
//
////    semi_local::hybrid_approach::hybrid(a,a_size,b,b_size,map,*actual,80000*2,3,2);
////       auto     actual =   semi_local::steady_ant_approach::get_semi_local_kernel(a,a_size,b,b_size,map);
//
//
//    auto time3 = std::chrono::high_resolution_clock::now() - begin3;
//    std::cout << "stikcy braid with steady ant approach " << std::chrono::duration<double, std::milli>(time3).count()
//              << std::endl;
//
//    actual->print(std::cout);


//
//    for (int i = 0;  i < a_size+b_size;i++){
//        should.set_point(res[i],i);
//    }

//    auto a1 = std::vector<int>(seq_a.begin(),seq_a.begin()+a_size/2);
//    auto a2 = std::vector<int>(seq_a.begin()+a_size / 2,seq_a.end());
//    auto l1 = semi_local::strand_combing_approach::sticky_braid_mpi<int,int>();



//    std::cout<<should.is_equal_to(*actual);

//std::cout<<std::endl;

//    should.print(std::cout);
//    delete actual;


}


//int main(int argc, char *argv[]) {
//    using namespace boost::multiprecision;
//
//    typedef  int  semi_type;
//     typedef  semi_type  wordType;
////todo fix description about required size
//
//    std::srand(0); // use current time as seed for random generator
//
//////
//////    auto a3 = new uint256_t[5];
////    a3[0] = 4;
////
////    std::cout<<a3[0];
////    std::cout<<sizeof (unsigned short)<<std::endl;
//////    int a_size = (4)*2000*3;//
////    int b_size = (6)*2000*3;
//
//
//    int a_size = 10*245;//
//    int b_size = 256*32+31;
//
//
//    auto seq_a = gen_vector_seq(a_size ,8);
//    auto seq_b = gen_vector_seq(b_size,8);
//
//    semi_type* a1 = new semi_type[a_size];
//    semi_type* b1 = new semi_type[b_size];
//
//    auto a2 = new int[a_size];
//    auto b2 = new int[b_size];
//
//    for (int i =0; i<a_size;i++) {
//        std::cout<<seq_a[i];
//        a1[i] = seq_a[i];
//        a2[i] = seq_a[i];
//    }
//    std::cout<<std::endl;
//    for (int i =0; i<b_size;i++) {
//        std::cout<<seq_b[i];
//        b1[i] = seq_b[i];
//        b2[i] = seq_b[i];
//    }
//    std::cout<<std::endl;
//    auto set = std::unordered_set<int>({0,1,2,3,4,5,6,7});
//    auto lookup = std::unordered_map<int,int *>();
//    lllcs_hyyro::preprocess(b2,b_size,set,lookup);
//
//    std::cout<<"siuzeof"<<sizeof (uint16_t)*8<<std::endl;
//
//    std::cout<<"hyyro:"<<lllcs_hyyro::hyyro_magic_mpi_with_precalc(a2,a_size,b_size/32+1,lookup,8)<<std::endl;
//    std::cout<<"hyyro_corre:"<<lllcs_hyyro::hyyro_magic(a2,a_size,b_size/(32 )+1,lookup)<<std::endl;
//
//    auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({7,6,5,4,3,2,1,0}));
//    auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second,32);
//    auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second,32);
//
//    for (int i =0; i<a.first.second;i++) {
//        std::cout<<std::bitset<8>(a.first.first[i])<<',';
//    }
//    std::cout<<std::endl;
//    for (int i =0; i<b.first.second;i++) {
//        std::cout<<std::bitset<8>(b.first.first[i])<<',';
//    }
//    std::cout<<std::endl;
//
//
//
////    mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({1,0}));
////    a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
////    b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
//
//
////
//
//
////
////    auto begin0 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << std::endl<<"res: "<< prefix_lcs_via_braid_bits_2symbol_v3_full_mask<unsigned int,wordType>
////            (a.first.first,a.first.second, a.second ,b.first.first,b.first.second, b.second,1) << std::endl;
////    auto time0 = std::chrono::high_resolution_clock::now() - begin0;
////    std::cout <<"time 2: " <<std::chrono::duration<double, std::milli>(time0).count() << std::endl;
//
//
////    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << std::endl<<"res: "<<prefix_lcs_via_semi_local::binary::llcs_2symbol_naive_combing(
////            a.first.first,a.first.second,b.first.first,b.first.second,a.second,1)<<std::endl;
////    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
////    std::cout <<"naive: " <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;
////
////    begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << std::endl<<"res: "<<prefix_lcs_via_semi_local::binary::llcs_2symbol_smart_combing(
////            a.first.first,a.first.second,b.first.first,b.first.second,a.second,1, true)<<std::endl;
////    time1 = std::chrono::high_resolution_clock::now() - begin1;
////    std::cout <<"true " <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;
////
////
//
////    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << std::endl<<"res: "<<prefix_lcs_via_semi_local::binary::llcs_2symbol_smart_combing(
////            a.first.first,a.first.second,b.first.first,b.first.second,a.second,1)<<std::endl;
////    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
////    std::cout <<"false" <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;
////
////    begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << std::endl<<"res: "<<prefix_lcs_via_semi_local::nary::llcs_nary_symbol_smart_combing(
////            a.first.first,a.first.second,b.first.first,b.first.second,a.second,5,1)<<std::endl;
////    time1 = std::chrono::high_resolution_clock::now() - begin1;
////    std::cout <<"nary:" <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//
////    auto perm = Permutation(a_size+b_size,a_size+b_size);
////    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    semi_local::strand_combing_approach::sticky_braid_mpi_branchless(perm,a1,a_size,b1,b_size);
////    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
////    std::cout <<"col="<< perm.get_col_by_row(0) <<"\n"<<"type:" <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//
//
//
//    auto begin4 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"res corret: "<< prefix_lcs_sequential(a1,a_size,b1,b_size) << std::endl;
//    auto time4 = std::chrono::high_resolution_clock::now() - begin4;
//    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time4).count() << std::endl;
//
//
//
//    //
////
////    std::cout << std::endl<<"Res binary mpi : "<< prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
//////                                                                                   b.first.first,b.first.second, b.second,2,thds) << std::endl;
////
////
////
////
////
////auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//
////    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
////    std::cout <<"Time prefix_lcs:" <<std::chrono::duration<double, std::milli>(time2).count() << std::endl;
////
////
////
////
////
////
//////
//////
////    auto begin4 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << std::endl<<"ResWithoutIF: "<< sticky_braid_sequential_without_if<int,int>(seq_a,seq_b) << std::endl;
////    auto time4 = std::chrono::high_resolution_clock::now() - begin4;
////    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time4).count() << std::endl;
////
////    auto begin5 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    std::cout << std::endl<<"ResWithif: "<< sticky_braid_sequential<int,int>(seq_a,seq_b) << std::endl;
////    auto time5 = std::chrono::high_resolution_clock::now() - begin5;
////    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time5).count() << std::endl;
//
//////
//////    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////////    prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second ,
////////                                         b.first.first,b.first.second, b.second,2,thds)
//////    std::cout <<"Res:"<< prefix_lcs_sequential(seq_a,seq_b) << std::endl;
//////    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//////    std::cout << std::chrono::duration<double, std::milli>(time2).count() << std::endl;
//////
//////        auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//////    std::cout << sticky_braid_mpi(seq_a, seq_b,1) << std::endl;
//////    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//////    std::cout << std::chrono::duration<double, std::milli>(time1).count() << std::endl;
////////
//////    auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//////    std::cout << sticky_braid_sequential (seq_a, seq_b) << std::endl;
//////    auto time3 = std::chrono::high_resolution_clock::now() - begin3;
//////    std::cout << std::chrono::duration<double, std::milli>(time3).count() << std::endl;
//////
////
////
////    //       ---
////    //     //   \\
////    //    //     \\
////    //   //=======\\
////    //  //         \\
////    //              ||
////    //
////
////
//    return 0 ;
//
//}
//
//
//
