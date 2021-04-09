
#include <string>
#include <iostream>
#include <chrono>
#include "../parsers.h"
#include "../naive_prefix_lcs.h"
#include "../prefix_lcs/bitwise/encoders_and_decoders.h"
#include "../prefix_lcs/bitwise/transposition_network_binary_alphabet.h"

/**
 * Accepts only binary strings --- 0 and 1 are allowed. Both should be divisible by 32
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
    int * _a = split(name_content_a.second.second,",",a_size);
    int * _b = split(name_content_b.second.second,",",b_size);

    auto a_vector = new std::vector<int>();
    auto b_vector = new std::vector<int>();
    for (int i = 0; i < a_size ; ++i) a_vector->push_back(_a[i]);
    for (int i = 0; i < b_size ; ++i) b_vector->push_back(_b[i]);

    auto mappers = encode_alphabet<int,unsigned int>(std::unordered_set<int>({0,1}));
    auto a = encode_reverse<int, unsigned int>(*a_vector, &mappers.first, &mappers.second,2);
    auto b = encode<int, unsigned int>(*b_vector,&mappers.first,&mappers.second,2);

    auto beg = std::chrono::high_resolution_clock::now();
    auto score =  prefix_lcs_via_semi_local::binary::llcs_2symbol_smart_combing<unsigned int>(
            a.first.first,a.first.second, b.first.first,b.first.second, a.second, thds,false);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());
    std::cout << 0 <<  std::endl; // some preprocess
    std::cout << elapsed_time << std::endl; // algo time
    std::cout << score << std::endl;
    std::cout<< a_size<<std::endl;
    std::cout<< b_size<<std::endl;
    std::cout<< a_name<<std::endl;
    std::cout<< b_name;

    delete a_vector;
    delete b_vector;
    delete _a;
    delete _b;

}

