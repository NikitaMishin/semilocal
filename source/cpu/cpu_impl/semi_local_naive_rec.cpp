//
////
//// Created by garrancha on 23.01.2021.
////
//
//
//
//#include <string>
//#include <iostream>
//#include <chrono>
//#include "../semi_local.h"
//#include "../parsers.h"
//#include "../unit_monge_mult/steady_ant.h"
//
//
///**
// * Solves semi-local problem for strings a and b.
// * Parallel  combing version with combined 1 and third phases.
// *  Use simd parallelism with antidiagonal pattern and thread-level parallelism
// * @param argc
// * @param argv
// * @return
// */
//int main(int argc, char *argv[]) {
//    int thds = strtol(argv[1], NULL, 10);
//    std::string a_filepath = std::string(argv[2]);
//    std::string b_filepath = std::string(argv[3]);
//    auto name_content_a = parse_input_file(a_filepath);
//    auto name_content_b = parse_input_file(b_filepath);
//
//    int a_size = name_content_a.first;
//    int b_size = name_content_b.first;
//    auto a_name = name_content_a.second.first;
//    auto b_name = name_content_b.second.first;
//    int * a = split(name_content_a.second.second,",",a_size);
//    int * b = split(name_content_b.second.second,",",b_size);
//
//
//
//    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
//    auto beg_precalc = std::chrono::high_resolution_clock::now();
//    precalc(map,5);
//    auto delta = std::chrono::high_resolution_clock::now() - beg_precalc;
//    auto precalc_elapsed_time = long(std::chrono::duration<double, std::milli>(delta).count());
//
//
//
//
//    auto beg = std::chrono::high_resolution_clock::now();
//    auto perm = semi_local::steady_ant_approach::get_semi_local_kernel(a,a_size,b,b_size,map);
//    auto time = std::chrono::high_resolution_clock::now() - beg;
//    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());
//    std::cout << precalc_elapsed_time <<   std::endl; // some preprocess
//    std::cout << elapsed_time <<  std::endl; // algo time
//    std::cout << hash(*perm, perm->row_size) << std::endl;
//    std::cout<< a_size<<std::endl;
//    std::cout<< b_size<<std::endl;
//    std::cout<< a_name<<std::endl;
//    std::cout<< b_name;
//
//    delete[] a;
//    delete[] b;
//}
////3336066023
