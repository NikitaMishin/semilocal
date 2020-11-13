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

void
precalc(std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>> &map,
        int max_size) {
    using namespace std;
    int p_arr[max_size];
    int q_arr[max_size];

    // f(p,q) = result//

    //   n x n matrices
    // calculcate [n/4 x n/4] parallel // O(n^2/4)
    //
    // log

    for (int size = 1; size < max_size + 1; size++) {


        for (int i = 0; i < size; ++i) p_arr[i] = i;


        do {
            for (int i = 0; i < size; ++i) q_arr[i] = i;

            do {
                auto p = new Permutation(size, size);
                for (int i = 0; i < size; ++i) p->set_point(i, p_arr[i]);

                auto q = new Permutation(size, size);
                for (int i = 0; i < size; ++i) q->set_point(i, q_arr[i]);
                long long hash_p = hash<AbstractPermutation>()(*p);
                long long hash_q = hash<AbstractPermutation>()(*q);

                auto r = steady_ant::steady_ant(p, q);
                auto points = std::vector<std::pair<int, int>>();
                if (map[size][hash_p].count(hash_q) > 0) {
                    std::cout << " Some error";
                    return;
                }
                r->to_points_on_grid(points);
                map[size][hash_p][hash_q] = points;

            } while (std::next_permutation(q_arr, q_arr + size));


        } while (std::next_permutation(p_arr, p_arr + size));
    }
}

int main() {
    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();

    auto a_size = 100;
    auto b_size = 100;

    auto seq_a = gen_vector_seq(a_size,25);
    auto seq_b = gen_vector_seq(b_size,25);
    int a[a_size];
    int b[b_size];
    for (int i = 0; i < a_size; ++i) {
        a[i] = seq_a[i];
    }
    for (int i = 0; i < b_size; ++i) {
        b[i] = seq_b[i];
    }

    //
    // |\ /|



    std::cout << "Precalc for value 5" << std::endl;
    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    precalc(map, 5);
    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
    std::cout << "Precalc " << std::chrono::duration<double, std::milli>(time1).count() << std::endl;

    std::cout << "Started on |a|="<<a_size <<", |b|="<<b_size << std::endl;

    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto res = sticky_braid_sequential<int,int>(seq_a,seq_b);
    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
    std::cout << "sticky braid with reducing approach " << std::chrono::duration<double, std::milli>(time2).count()
              << std::endl;

    auto should = Permutation(a_size+b_size,a_size+b_size);
    for (int i = 0;  i < a_size+b_size;i++){
        should.set_point(res[i],i);
    }

    auto a1 = std::vector<int>(seq_a.begin(),seq_a.begin()+a_size/2);
    auto a2 = std::vector<int>(seq_a.begin()+a_size / 2,seq_a.end());
    auto l1 = sticky_braid_sequential<int,int>(a1,seq_b);




    auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false


    auto actual = semi_local_lcs::get_semi_local_kernel(a,a_size,b,b_size,map);


    auto time3 = std::chrono::high_resolution_clock::now() - begin3;
    std::cout << "stikcy braid with steady ant approach " << std::chrono::duration<double, std::milli>(time3).count()
              << std::endl;


    std::cout<<should.is_equal_to(*actual);




}

AbstractPermutation *run(AbstractPermutation *m, AbstractPermutation *n, int *memory_block, PrecalcMap &map) {
    auto m_new = new PermutationPreAllocated(n->row_size, n->col_size, memory_block, memory_block + n->row_size);
    auto n_new = PermutationPreAllocated(n->row_size, n->col_size, memory_block + 2 * n->row_size,
                                         memory_block + 3 * n->row_size);

    for (int i = 0; i < m->row_size; ++i) {
        auto col = m->get_col_by_row(i);
        m_new->set_point(i, col);
    }

    for (int i = 0; i < n->row_size; ++i) {
        auto col = n->get_col_by_row(i);
        n_new.set_point(i, col);
    }

    int nearest_2_degree = pow(2, int(ceil(log2(2 * n->row_size))));
    int total = int(log2(nearest_2_degree)) * nearest_2_degree;


    //8n memory
    steady_ant::steady_ant_with_precalc_and_memory(m_new, &n_new, memory_block,
                                                   memory_block + 4 * n->row_size,
                                                   memory_block + 8 * n->row_size, map, total);


    return m_new;


}

//int main(int argc, char *argv[]) {
//    omp_set_max_active_levels(2);
//    std::cout<< omp_get_active_level();
//    auto m = distance_product::get_permutation_matrix(10, 10, -2);
////    m->print(std::cout);
////10mil
//    auto n = distance_product::get_permutation_matrix(10, 10, -3);
//
//
//
//    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
//
////    auto expected = distance_product::naive::mult_dist(m, n);
////
////    std::cout<<sizeof(n);
//
////    omp_set_nested(0);
//    omp_set_num_threads(4);
//    std::cout << "Precalc for value 5" << std::endl;
//    auto begin1 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    precalc(map, 5);
//    auto time1 = std::chrono::high_resolution_clock::now() - begin1;
//    std::cout << "Precalc " << std::chrono::duration<double, std::milli>(time1).count() << std::endl;
//
//    std::cout << "Memory alloc for 20n" << std::endl;
//    auto begin4 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    int nearest_2_degree = pow(2,int(ceil(log2(2*n->row_size))));
//    auto memory_block = new int[n->row_size * 8 + int(log2(nearest_2_degree)) * nearest_2_degree];
//    auto time4 = std::chrono::high_resolution_clock::now() - begin4;
//    std::cout << "Time " << std::chrono::duration<double, std::milli>(time4).count() << std::endl;
/////row_to_col
//// col_to_row
////n/2, n/2
//// n/p
//// 2 ->n
//
//    auto product = Permutation(11,11);
//    steady_ant::staggered_sticky_multiplication(m,n,9,map,&product);
//    m->print(std::cout);
//    std::cout<<std::endl;
//    n->print(std::cout);
//    std::cout<<std::endl;
//    product.print(std::cout);
//    std::cout<<std::endl;
//
//
//    std::cout << "Started on 10000000x10000000" << std::endl;
//    auto begin2 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    auto actual = run(m, n, memory_block, map);
//    std::cout << std::endl << "res: " << actual->get_col_by_row(7) << std::endl;
//    auto time2 = std::chrono::high_resolution_clock::now() - begin2;
//    std::cout << "with precalc and memmory alloc " << std::chrono::duration<double, std::milli>(time2).count()
//              << std::endl;
//
////    m = distance_product::get_permutation_matrix(10000000, 10000000, -233);
//////    m->print(std::cout);
////
////    n = distance_product::get_permutation_matrix(10000000, 10000000, -2);
////
////
//    auto begin0 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
//    auto exp = steady_ant::steady_ant(m, n);
////    std::cout << std::endl<<"res: "<< steady_ant::steady_ant (m,n)->get_col_by_row(34) << std::endl;
//    auto time0 = std::chrono::high_resolution_clock::now() - begin0;
//    std::cout << "without precalc and without memory alloc: "
//              << std::chrono::duration<double, std::milli>(time0).count() << std::endl;
//
////    delete res;
//
//    std::cout << exp->is_equal_to(product);
//    delete[] memory_block;
//    delete n;
//    delete m;
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
