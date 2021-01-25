
//
// Created by garrancha on 23.01.2021.
//



#include <string>
#include <iostream>
#include <chrono>
#include "../semi_local.h"
#include "../fasta_parser.h"

/**
 * Solves semi-local problem for strings a and b.
 * Parallel  combing version with reverse storage of a.
 * Use simd parallelism with antidiagonal pattern and thread-level parallism
 * @param argc
 * @param argv
 * @return
 */
int main(int argc, char *argv[]) {
    int thds = strtol(argv[1], NULL, 10);
    std::string a_filepath = std::string(argv[2]);
    std::string b_filepath = std::string(argv[3]);
    auto name_content_a = parse_input_file(a_filepath);
    auto name_content_b = parse_input_file(b_filepath);

    int a_size = name_content_a.first;
    int b_size = name_content_b.first;
    auto a_name = name_content_a.second.first;
    auto b_name = name_content_b.second.first;
    int * a = split(name_content_a.second.second,",",a_size);
    int * b = split(name_content_b.second.second,",",b_size);


    auto beg_precalc = std::chrono::high_resolution_clock::now();


    int * a_reverse = new int[a_size];
    for (int i = 0; i < a_size; ++i) {
        a_reverse[i] = a[a_size - 1 - i];
    }

    auto delta = std::chrono::high_resolution_clock::now() - beg_precalc;
    auto precalc_elapsed_time = long(std::chrono::duration<double, std::milli>(delta).count());


    auto perm = Permutation(a_size+b_size,a_size+b_size);

    auto beg = std::chrono::high_resolution_clock::now();
    semi_local::strand_combing_approach::sticky_braid_mpi_reversea(perm, a_reverse, a_size, b, b_size, thds);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());
    std::cout << precalc_elapsed_time <<  "ms"  << std::endl; // some preprocess
    std::cout << elapsed_time << "ms" << std::endl; // algo time
    std::cout << hash(perm, perm.row_size) << std::endl;
    std::cout<< a_name<<std::endl;
    std::cout<< b_name<<std::endl;

    delete[] a;
    delete[] b;
    delete[] a_reverse;
}