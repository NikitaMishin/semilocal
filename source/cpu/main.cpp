//
// Created by nikita on 27.07.2020.
//
//
#include "semi_local.h"
#include "transposition_network_approach/encoders_and_decoders.h"
#include "sequence_generators.h"
#include "fasta_parser.h"
#include "naive_prefix_lcs.h"
#include "transposition_network_approach/4symbol_new.h"
#include "transposition_network_approach/2symbol_new.h"
#include "transposition_network_approach/2symbol_new_2.h"
#include "unit_monge_mult/steady_ant.h"
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



int main() {

    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
    omp_set_nested(true);
    auto a_size = 25000/2;
    auto b_size = 25000*4;

    auto seq_a = gen_vector_seq(a_size,40);
    auto seq_b = gen_vector_seq(b_size,40);

    int a[a_size];
    int b[b_size];
    for (int i = 0; i < a_size; ++i) {
        a[i] = seq_a[i];
//        std::cout<<a[i]<<" ";
    }
    for (int i = 0; i < b_size; ++i) {
        b[i] = seq_b[i];
//        std::cout<<b[i]<<" ";
    }
    std::cout<<std::endl;

    //
    // |\ /|



    std::cout << "Precalc for value 5" << std::endl;
    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    distance_unit_monge_product::steady_ant::precalc(map, 5);
    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
    std::cout << "Precalc " << std::chrono::duration<double, std::milli>(time1).count() << std::endl;

    std::cout << "Started on |a|="<<a_size <<", |b|="<<b_size << std::endl;





    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false

    auto should = Permutation(a_size+b_size,a_size+b_size);




    semi_local::hybrid_approach::first_and_third_phase_merged(a, a_size, b, b_size, should, map, 2,32);
    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
    std::cout << "merged " << std::chrono::duration<double, std::milli>(time2).count()
              << std::endl;


    auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false

    auto actual = new  Permutation(a_size+b_size,a_size+b_size);
    semi_local::strand_combing_approach::sticky_braid_mpi(*actual,a,a_size,b,b_size,32);

//    semi_local::hybrid_approach::hybrid(a,a_size,b,b_size,map,*actual,80000*2,3,2);
//       auto     actual =   semi_local::steady_ant_approach::get_semi_local_kernel(a,a_size,b,b_size,map);


    auto time3 = std::chrono::high_resolution_clock::now() - begin3;
    std::cout << "stikcy braid with steady ant approach " << std::chrono::duration<double, std::milli>(time3).count()
              << std::endl;


//
//    for (int i = 0;  i < a_size+b_size;i++){
//        should.set_point(res[i],i);
//    }

//    auto a1 = std::vector<int>(seq_a.begin(),seq_a.begin()+a_size/2);
//    auto a2 = std::vector<int>(seq_a.begin()+a_size / 2,seq_a.end());
//    auto l1 = semi_local::strand_combing_approach::sticky_braid_mpi<int,int>();



    std::cout<<should.is_equal_to(*actual);

//std::cout<<std::endl;

//    should.print(std::cout);
//    delete actual;


}


//int main(int argc, char *argv[]) {
//    omp_set_nested(true);
//    omp_set_dynamic(true);
////    omp_set_num_threads(16);
////    omp_set
////omp_set_max_active_levels(2);
//    std::cout<< omp_get_num_threads();
//    auto size = 2500000;
//
//    auto m = new Permutation(size,size);
//    fill_permutation_matrix(m,size, size, -2);
//
//    auto n = new Permutation(size,size);
//    fill_permutation_matrix(n,size, size, -3);
//
//
//
////res.print(std::cout);
//    std::cout<<std::endl;
////    exp->print(std::cout);
////    delete exp;
////    auto expected = distance_product::naive::mult_dist(m, n);
////
////    std::cout<<sizeof(n);
//
////    omp_set_num_threads(4);
//    std::cout << "Precalc for value 5" << std::endl;
//    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
//    distance_unit_monge_product::steady_ant::precalc(map, 5);
//    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//    std::cout << "Precalc " << std::chrono::duration<double, std::milli>(time1).count() << std::endl;
////
//    std::cout << "Memory alloc for 20n" << std::endl;
//    auto begin4 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    int nearest_2_degree = pow(2,int(ceil(log2(2*n->row_size))));
//    auto memory_block = new int[n->row_size * 8 + int(log2(nearest_2_degree)) * nearest_2_degree];
//    auto time4 = std::chrono::high_resolution_clock::now() - begin4;
//    std::cout << "Time " << std::chrono::duration<double, std::milli>(time4).count() << std::endl;
//
//    std::cout << "Started on 10000000x10000000" << std::endl;
//    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    auto actual = run(m, n, memory_block, map);
//    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    std::cout << "with precalc and memmory alloc " << std::chrono::duration<double, std::milli>(time2).count()
//              << std::endl;
////
//////    m = distance_product::get_permutation_matrix(10000000, 10000000, -233);
////////    m->print(std::cout);
//////
//////    n = distance_product::get_permutation_matrix(10000000, 10000000, -2);
//////
//////
////    auto begin0 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
////    auto exp = steady_ant::steady_ant(m, n);
//////    std::cout << std::endl<<"res: "<< steady_ant::steady_ant (m,n)->get_col_by_row(34) << std::endl;
////    auto time0 = std::chrono::high_resolution_clock::now() - begin0;
////    std::cout << "without precalc and without memory alloc: "
////              << std::chrono::duration<double, std::milli>(time0).count() << std::endl;
////
//////    delete res;
//
//    std::cout << "Started on 10000000x10000000" << std::endl;
//    auto begin5 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    auto exp = distance_unit_monge_product::steady_ant::steady_ant(m,n);
//
//    auto time5 = std::chrono::high_resolution_clock::now() - begin5;
//    std::cout << "withOUT ANY " << std::chrono::duration<double, std::milli>(time5).count()
//              << std::endl;
//
//    std::cout << exp->is_equal_to(*actual);
////    delete[] memory_block;
////    delete n;
////    delete m;
//}


//
//int main(int argc, char *argv[]) {
//
//    typedef  unsigned int wordType;
//
//
////    std::cout<<sizeof(wordType) <<std::endl;
//    std::srand(std::time(NULL)); // use current time as seed for random generator
//
//
//
//
//
//    int a_size = 64*1000;//
//    int b_size = 64*1000;
////    int a_size = 1024*32*10;//
////    int b_size = 1024*1024*64;
////96458
//
//    auto seq_a = gen_vector_seq(a_size ,1);
//    auto seq_b = gen_vector_seq(b_size,1);
//
//
//    auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({3,2,1,0}));
//    auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
//    auto b = encode<int, wordType >(seq_b,&mappers.first,&mappers.second);
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
//    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"res: "<<
//              prefix_lcs_via_braid_bits_4symbol_v2_full_mask<wordType>(a.first.first,a.first.second, a.second ,
//                                                                       b.first.first,b.first.second, b.second,1)<< std::endl;
//    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//    std::cout <<"Time 4: " <<std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//
////    std::cout << std::endl<<"res corret: "<< sticky_braid_sequential_without_if<int,int>(seq_a,seq_b) << std::endl;
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
//    auto begin4 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"ResWithoutIF: "<< sticky_braid_sequential_without_if<int,int>(seq_a,seq_b) << std::endl;
//    auto time4 = std::chrono::high_resolution_clock::now() - begin4;
//    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time4).count() << std::endl;
//
//    auto begin5 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    std::cout << std::endl<<"ResWithif: "<< sticky_braid_sequential<int,int>(seq_a,seq_b) << std::endl;
//    auto time5 = std::chrono::high_resolution_clock::now() - begin5;
//    std::cout <<"Time:" <<std::chrono::duration<double, std::milli>(time5).count() << std::endl;
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
////    return 0 ;
//
//}
//
//
//
