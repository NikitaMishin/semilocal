//
// Created by nikita on 30.07.2020.
//

#ifndef CPU_UTILS_H
#define CPU_UTILS_H


#include <cmath>

/**
 *
 * @tparam Input
 * @tparam Output
 * @param a
 * @param n
 * @param total_symbols
 * @param encoder
 * @return
 */
template<class Input, class Output>
std::vector<Output> decode(Input const *a, int n, int total_symbols, std::unordered_map<Input, Output> encoder) {
    auto alphabet_size = encoder.size();
    auto bits_per_symbol = int(std::ceil(log2(alphabet_size)));
    auto shift = bits_per_symbol;
    auto word_size_in_bits = sizeof(Input) * 8;
    auto symbols_in_word = word_size_in_bits / shift;
    auto sequence = new std::vector<Output>();
    auto mask = (Input(1) << bits_per_symbol) - 1;

//    fill bitset
    for (int i = 0; i < n - 1; ++i) {
        auto word = a[i];
        for (int symbol = 0; symbol < symbols_in_word; symbol++) {
            sequence->push_back(encoder[(word >> shift * symbol) & mask]);
        }
    }

//    fill last
    for (int i = n - 1; i < n; ++i) {
        auto word = a[i];
        for (int symbol = 0; i * symbols_in_word + symbol < total_symbols; symbol++) {
            sequence->push_back(encoder[(word >> shift * symbol) & mask]);
        }
    }


//    for (auto &it: *mapper_reverse) {
//        std::cout << std::bitset<2>(it.first) << ":" << it.second << std::endl;
//    }
//    for (auto &it:a) {
//        std::cout << it;
//    }
//    std::cout << std::endl;

    return *sequence;


}


/**
 * Generate sequence of int numbers
 * @param size int
 * @param alphabet_size int, default = 4
 * @return vector<int> of size @size
 */
std::vector<int> gen_vector_seq(int size, int alphabet_size = 4) {

    auto v = std::vector<int>();
    for (int i = 0; i < size; ++i) {
        v.push_back(rand() % alphabet_size);
    }
    return v;
}

/**
 * Generate sequence of chars
 * @param size int
 * @return vector<int> of size @size
 */
std::vector<char> gen_vector_ACGT(int size) {

    auto v = std::vector<char>();
    for (int i = 0; i < size; ++i) {
        switch(int(rand()) % 4) {
            case 0:
                v.push_back('A');
                break;
            case 1:
                v.push_back('C');
                break;
            case 2:
                v.push_back('G');
                break;
            case 3:
                v.push_back('T');
                break;
        }
    }
    return v;
}


/**
 * Generate sequence of int numbers
 * @param size int
 * @param alphabet_size int, default = 4
 * @return array<int> of size @size
 */
 int * gen_array_seq(int size, int alphabet_size = 4) {

    auto v = new int[size];
    for (int i = 0; i < size; ++i) {
        v[i]  = rand() % alphabet_size;
    }
    return v;
}







#endif //CPU_UTILS_H
