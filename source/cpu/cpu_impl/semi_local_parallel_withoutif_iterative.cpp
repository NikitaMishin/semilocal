
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
 * Parallel  combing version without branching.
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

    auto perm = Permutation(a_size+b_size,a_size+b_size);

    auto beg = std::chrono::high_resolution_clock::now();
    semi_local::strand_combing_approach::sticky_braid_mpi_withoutif(perm, a, a_size, b, b_size, thds);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());
    std::cout << 0   << std::endl; // some preprocess
    std::cout << elapsed_time  << std::endl; // algo time
    std::cout << hash(perm, perm.row_size) << std::endl;
    std::cout<< a_size<<std::endl;
    std::cout<< b_size<<std::endl;
    std::cout<< a_name<<std::endl;
    std::cout<< b_name;


    delete[] a;
    delete[] b;

}