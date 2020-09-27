
#include <string>
#include <iostream>
#include <chrono>
#include "../fasta_parser.h"
#include "../naive_prefix_lcs.h"

//static const int length = 1024*8;

int main(int argc, char *argv[]) {
    int thds = strtol(argv[1], NULL, 10);
    std::string a_filepath = std::string(argv[2]);
    std::string b_filepath = std::string(argv[3]);
    auto name_content_a = parse_input_file(a_filepath);
    auto name_content_b = parse_input_file(b_filepath);

    auto seq_a = transform_to_int_vector(name_content_a.second.second);
    auto seq_b = transform_to_int_vector(name_content_b.second.second);
    auto beg = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
    auto res = prefix_lcs_sequential(seq_a, seq_b);
    auto time = std::chrono::high_resolution_clock::now() - beg;
    auto elapsed_time = long(std::chrono::duration<double, std::milli>(time).count());

    std::cout << 0 << std::endl;
    std::cout << elapsed_time << std::endl;
    std::cout << res << std::endl;
    std::cout << seq_a.size() << std::endl;
    std::cout << seq_b.size();

}
