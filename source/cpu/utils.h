//
// Created by nikita on 30.07.2020.
//

#ifndef CPU_UTILS_H
#define CPU_UTILS_H


#include <cmath>
#include <algorithm>
#include <unordered_set>


/**
 * Maps alphabet symbols to numbers from interval 0 to alpabet_size
 * @tparam Input  any comparable type and hashable
 * @tparam Output number type. Unsigned int/U long,size_t,short
 * @param alphabet_set
 * @return
 */
template<class Input, class Output>
std::pair<std::unordered_map<Input, Output>, std::unordered_map<Output, Input>>
encode_alphabet(std::unordered_set<Input> alphabet_set) {
    auto mapper_forward = new std::unordered_map<Input, Output>;
    auto mapper_reverse = new std::unordered_map<Output, Input>;
    auto encoder = Output(0);
    for (Input s : alphabet_set) {
        if (mapper_forward->count(s) == 0) {
            (*mapper_forward)[s] = Output(encoder);
            (*mapper_reverse)[encoder] = s;
            encoder++;
        }
    }
    return std::pair<std::unordered_map<Input, Output>, std::unordered_map<Output, Input>>(*mapper_forward,
                                                                                           *mapper_reverse);
}


/**
 * Encode given sequence with Input type symbols to bit array packed in Output type
 * according to alphabet_mapping and amount of bits per symbol
 * @tparam Input
 * @tparam Output
 * @param a
 * @return (packed array, amount of packed elements, size of initial sequence)
 */
template<class Input, class Output>
std::pair<std::pair<Output *, int>, int>
encode(std::vector<Input> const &a, std::unordered_map<Input, Output> *mapper_forward,
       std::unordered_map<Output, Input> *mapper_reverse) {

    auto alphabet_size = mapper_reverse->size();
    auto bits_per_symbol = int(std::ceil(log2(alphabet_size)));
    auto shift = bits_per_symbol;
    auto word_size_in_bits = sizeof(Output) * 8;
    auto symbols_in_word = int(word_size_in_bits / shift);

    auto bytes_needed = int(std::ceil(a.size() * 1.0 / symbols_in_word) * sizeof(Output));

    auto bitset_array = static_cast<Output *> (aligned_alloc(sizeof(Output), bytes_needed));
    auto n = bytes_needed / sizeof(Output);


//    fill bitset
    for (int i = 0; i < n - 1; ++i) {
        Output word = 0;
        for (int symbol = 0; symbol < symbols_in_word; symbol++) {
            word |= ((*mapper_forward)[a[i * symbols_in_word + symbol]]) << shift * symbol;
        }
        bitset_array[i] = word;
    }

//    fill last
    for (int i = n - 1; i < n; ++i) {
        Output word = 0;
        for (int symbol = 0; i * symbols_in_word + symbol < a.size(); symbol++) {
            word |= (*mapper_forward)[a[i * symbols_in_word + symbol]] << shift * symbol;
        }
        bitset_array[i] = word;
    }

    return std::make_pair(std::make_pair(bitset_array, n), a.size());
}

/*
*
* Encode given sequence with Input type symbols to  reversed bit array packed in Output type
* according to alphabet_mapping and amount of bits per symbol.
 * I.e abcde -> edcba -> to bits for each symbol
 * @tparam Input
 * @tparam Output
 * @param a
 * @return (packed array, amount of packed elements, size of initial sequence)
 */
template<class Input, class Output>
std::pair<std::pair<Output *, int>, int>
encode_reverse(std::vector<Input> const &a, std::unordered_map<Input, Output> *mapper_forward,
               std::unordered_map<Output, Input> *mapper_reverse) {
    return encode<Input, Output>(std::vector<Input>(a.rbegin(), a.rend()), mapper_forward, mapper_reverse);
}


/**
 * Decode given packed sequence in bits in numeric type Input to symbol vector of type Output
 * according to decoder and amount of bits per symbol
 * @tparam Input
 * @tparam Output
 * @param a
 * @param n
 * @param total_symbols
 * @param decoder
 * @return
 */
template<class Input, class Output>
std::vector<Output> decode(Input const *a, int n, int total_symbols, std::unordered_map<Input, Output> decoder) {
    auto alphabet_size = decoder.size();
    auto bits_per_symbol = int(std::ceil(log2(alphabet_size)));
    auto shift = bits_per_symbol;
    auto word_size_in_bits = sizeof(Input) * 8;
    auto symbols_in_word = int(word_size_in_bits / shift);
    auto sequence = new std::vector<Output>();
    auto mask = (Input(1) << bits_per_symbol) - 1;

//    fill bitset
    for (int i = 0; i < n - 1; ++i) {
        auto word = a[i];
        for (int symbol = 0; symbol < symbols_in_word; symbol++) {
            sequence->push_back(decoder[(word >> shift * symbol) & mask]);
        }
    }

//    fill last
    for (int i = n - 1; i < n; ++i) {
        auto word = a[i];
        for (int symbol = 0; i * symbols_in_word + symbol < total_symbols; symbol++) {
            sequence->push_back(decoder[(word >> shift * symbol) & mask]);
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
 * Decode given reversed packed sequence in bits in numeric type Input to symbol vector of type Output
 * according to decoder and amount of bits per symbol
 * @tparam Input
 * @tparam Output
 * @param a
 * @param n
 * @param total_symbols
 * @param decoder
 * @return
 */
template<class Input, class Output>
std::vector<Output>
decode_reverse(Input const *a, int n, int total_symbols, std::unordered_map<Input, Output> decoder) {
    auto v = decode(a, n, total_symbols, decoder);
    std::reverse(v.begin(), v.end());
    return v;
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
 * Generate sequence of chars of A C G T symbols
 * @param size int
 * @return vector<int> of size @size
 */
std::vector<char> gen_vector_ACGT(int size) {

    auto v = std::vector<char>();
    for (int i = 0; i < size; ++i) {
        switch (int (rand()) % 4) {
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
 * Generate sequence of int numbers with |alphabet|=4
 * @param size int
 * @param alphabet_size int, default = 4
 * @return array<int> of size @size
 */
int *gen_array_seq(int size, int alphabet_size = 4) {

    auto v = new int[size];
    for (int i = 0; i < size; ++i) {
        v[i] = rand() % alphabet_size;
    }
    return v;
}


#endif //CPU_UTILS_H
