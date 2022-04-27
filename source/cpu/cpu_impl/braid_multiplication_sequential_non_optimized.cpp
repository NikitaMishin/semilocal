
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
 * Program that evaluates performance of non-optimized sequential braid multiplication without precalc
 *
 * @param argc
 * @param argv Accepts  a dummy parameter depth , N - matrix size, and seed
 * @return std out in specified format
 */
int main(int argc, char *argv[]) {
    int depth = strtol(argv[1], NULL, 10);
    int n = strtol(argv[2], NULL, 10);
    int seed = strtol(argv[3], NULL, 10);


    auto p = new Permutation(n,n);
    auto q = new Permutation(n,n);

    fillPermutationMatrix(p, n, n, seed);
    fillPermutationMatrix(q, n, n, -seed);




    auto empty_map = PrecalcMap();

    auto beg = std::chrono::high_resolution_clock::now();
    auto product = distance_unit_monge_product::steady_ant::steady_ant(p,q,empty_map);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());


    std::cout << 0 <<  std::endl; // some preprocess
    std::cout << elapsed_time << std::endl; // algo time
    std::cout << hash(*product, product->row_size) << std::endl;
    std::cout<< n<<std::endl;
    std::cout<< seed;

}
