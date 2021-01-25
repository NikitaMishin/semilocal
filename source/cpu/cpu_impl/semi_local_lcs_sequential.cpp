
#include <string>
#include <iostream>
#include <chrono>
#include "../semi_local.h"
#include "../fasta_parser.h"

const long long R = 4294967279;
const long long M = 4294967291;

//static const int length = 1024*8;
template<class Input>
long long hash(Input *arr, int size) {
    long long hash = 0;
    for (int i = 0; i < size; i++)
        hash = (R * hash + arr[i]) % M;
    return hash;
}

int main(int argc, char *argv[]) {
    int thds = strtol(argv[1], NULL, 10);
    std::string a_filepath = std::string(argv[2]);
    std::string b_filepath = std::string(argv[3]);
    auto name_content_a = parse_input_file(a_filepath);
    auto name_content_b = parse_input_file(b_filepath);

    auto seq_a = transform_to_int_vector(name_content_a.second.second);
    auto seq_b = transform_to_int_vector(name_content_b.second.second);
    auto beg = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto res = semi_local::strand_combing_approach::sticky_braid_sequential_without_if<int, int>(seq_a, seq_b);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());

    std::cout << 0 << std::endl;
    std::cout << elapsed_time << std::endl;
    std::cout << hash<int>(res, seq_b.size() + seq_a.size()) << std::endl;
    std::cout << seq_a.size() << std::endl;
    std::cout << seq_b.size();
}