
//
// Created by garrancha on 23.01.2021.
//



#include <string>
#include <iostream>
#include <chrono>
#include "../semi_local.h"
#include "../fasta_parser.h"



int main(int argc, char *argv[]) {
    int depth = strtol(argv[1], NULL, 10);
    int n = strtol(argv[2], NULL, 10);
    int seed = strtol(argv[3], NULL, 10);


    auto p = new Permutation(n,n);
    auto q = new Permutation(n,n);

    fill_permutation_matrix(p,n,n,seed);
    fill_permutation_matrix(q,n,n,-seed);



    auto beg = std::chrono::high_resolution_clock::now();
    auto product = distance_unit_monge_product::steady_ant::steady_ant(p,q);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());


    std::cout << 0 <<  "ms"  << std::endl; // some preprocess
    std::cout << elapsed_time << "ms" << std::endl; // algo time
    std::cout << hash(*product, product->row_size) << std::endl;
    std::cout<< n<<std::endl;
    std::cout<< seed<<std::endl;

}
