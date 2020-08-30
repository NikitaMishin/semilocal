
#include <string>
#include <iostream>
#include <chrono>
#include "../fasta_parser.h"
#include "../transposition_network_approach/encoders_and_decoders.h"
#include "../transposition_network_approach/transposition_network_binary_alphabet.h"

//static const int length = 1024*8;

int main(int argc, char *argv[])
{
    typedef long long wordType;
    int thds = strtol(argv[1], NULL, 10);
    std::string a_filepath = std::string(argv[2]);
    std::string b_filepath = std::string(argv[3]);
    auto name_content_a = parse_input_file(a_filepath);
    auto name_content_b = parse_input_file(b_filepath);

    auto seq_a = transform_to_int_vector(name_content_a.second.second);
    auto seq_b = transform_to_int_vector(name_content_b.second.second);
    if(seq_a.size()>seq_b.size()) std::swap(seq_a,seq_b);

    auto beg_preprocess = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto mappers = encode_alphabet<int,wordType>(std::unordered_set<int>({1,0}));
    auto a = encode_reverse<int, wordType>(seq_a,&mappers.first,&mappers.second);
    auto b = encode<int, wordType>(seq_b,&mappers.first,&mappers.second);
    auto time_preprocess = std::chrono::high_resolution_clock::now() - beg_preprocess;
    auto elapsed_time_preprocess  = long (std::chrono::duration<double, std::milli>(time_preprocess).count());


    auto beg = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto res = prefix_lcs_via_braid_bits_binary_mpi(a.first.first,a.first.second, a.second, b.first.first,b.first.second, b.second, thds);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time  = long (std::chrono::duration<double, std::milli>(time).count());

    std::cout<<elapsed_time_preprocess<<std::endl;
    std::cout<<elapsed_time<<std::endl;
    std::cout<<res<<std::endl;
    std::cout<<seq_a.size()<<std::endl;
    std::cout<<seq_b.size();

}

