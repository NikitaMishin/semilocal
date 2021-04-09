
//
// Created by garrancha on 23.01.2021.
//



#include <string>
#include <iostream>
#include <chrono>
#include "../semi_local.h"
#include "../parsers.h"
#include "../test_utils.h"
//SET OPTIMAL THRESHOLD


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
    int *a = split(name_content_a.second.second, ",", a_size);
    int *b = split(name_content_b.second.second, ",", b_size);

    auto b_sp = new unsigned short[b_size];
    auto a_sp = new unsigned short[a_size];

    for (int j = 0; j < a_size; ++j) a_sp[j] = a[j];
    for (int j = 0; j < b_size; ++j) b_sp[j] = b[j];


    // get raz


    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
    auto beg_precalc = std::chrono::high_resolution_clock::now();
    precalc(map, 5);
    auto delta = std::chrono::high_resolution_clock::now() - beg_precalc;
    auto precalc_elapsed_time = long(std::chrono::duration<double, std::milli>(delta).count());


    auto perm = Permutation(a_size + b_size, a_size + b_size);

    auto beg = std::chrono::high_resolution_clock::now();
    semi_local::hybrid_iterative_wrapper<false>(perm, a_sp, a_size, b_sp, b_size, map, thds);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());
    std::cout << precalc_elapsed_time << std::endl; // some preprocess
    std::cout << elapsed_time << std::endl; // algo time
    std::cout << hash(perm, perm.row_size) << std::endl;
    std::cout << a_size << std::endl;
    std::cout << b_size << std::endl;
    std::cout << a_name << std::endl;
    std::cout << b_name;

    delete[] a;
    delete[] b;
    delete[] a_sp;
    delete[] b_sp;
}