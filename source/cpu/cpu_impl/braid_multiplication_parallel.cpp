
//
// Created by garrancha on 23.01.2021.
//



#include <string>
#include <iostream>
#include <chrono>
#include "../semi_local.h"
#include "../parsers.h"
#include "../test_utils.h"


/**
 * Program that evaluates performance of optimized parallel braid multiplication with precalc
 *
 * @param argc
 * @param argv Accepts depth that determines how many upper levels of recursion should be parallelized, N - matrix size,  and seed
 * @return std out in specified format
 */
int main(int argc, char *argv[]) {
    int depth = strtol(argv[1], NULL, 10);
    int n = strtol(argv[2], NULL, 10);
    int seed = strtol(argv[3], NULL, 10);


    omp_set_nested(true);

    auto p = Permutation(n,n);
    auto q = Permutation(n,n);
    auto product = Permutation(n,n);

    fill_permutation_matrix(&p,n,n,seed);
    fill_permutation_matrix(&q,n,n,-seed);

    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
    auto beg_precalc = std::chrono::high_resolution_clock::now();
    precalc(map,5);
    auto delta = std::chrono::high_resolution_clock::now() - beg_precalc;
    auto precalc_elapsed_time = long(std::chrono::duration<double, std::milli>(delta).count());




    auto beg = std::chrono::high_resolution_clock::now();
    distance_unit_monge_product::steady_ant::steady_ant_parallel_wrapper(p, q, product, map, depth);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());


    std::cout << precalc_elapsed_time << std::endl; // some preprocess
    std::cout << elapsed_time << std::endl; // algo time
    std::cout << hash(product, product.row_size) << std::endl;
    std::cout<< n<<std::endl;
    std::cout<< seed;

}
