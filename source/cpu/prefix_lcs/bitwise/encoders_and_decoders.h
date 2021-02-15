//
// Created by nikita on 19.08.2020.
//

#ifndef CPU_ENCODERS_AND_DECODERS_H
#define CPU_ENCODERS_AND_DECODERS_H

#include <cmath>
#include <algorithm>
#include <unordered_set>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <bitset>

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
    std::unordered_map<Input, Output> * mapper_forward = new std::unordered_map<Input, Output>;
    std::unordered_map<Output, Input> * mapper_reverse = new std::unordered_map<Output, Input>;
    Output encoder = Output(0);
    for (Input s : alphabet_set) {
        if (mapper_forward->count(s) == 0) {
            mapper_forward->insert(std::make_pair(s,encoder));
            mapper_reverse->insert(std::make_pair(encoder,s));
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


    auto alphabet_size = mapper_reverse->size();
    auto bits_per_symbol = int(std::ceil(log2(alphabet_size)));
    auto shift = bits_per_symbol;
    auto word_size_in_bits = sizeof(Output) * 8;
    auto symbols_in_word = int(word_size_in_bits / shift);

    auto bytes_needed = int(std::ceil(a.size() * 1.0 / symbols_in_word) * sizeof(Output));

    auto bitset_array = static_cast<Output *> (aligned_alloc(sizeof(Output), bytes_needed));
    auto n = bytes_needed / sizeof(Output);


    // fill first guy
    for (int i = 0; i < n-1 ; ++i) {
        Output word = Output(0);
        for (int symbol = 0; symbol < symbols_in_word; symbol++) {
            word |= ((*mapper_forward)[a[i * symbols_in_word + symbol]]) <<  (bits_per_symbol*(symbols_in_word - symbol - 1));
            }
        bitset_array[n-1-i] = word;
    }

    //    fill last
    for (int i = n - 1; i < n; ++i) {
        Output word = 0;
        for (int symbol = 0; (n-1) * symbols_in_word + symbol < a.size(); symbol++) {
            word |= ((*mapper_forward)[a[i * symbols_in_word + symbol]]) <<  (bits_per_symbol*(symbols_in_word - symbol - 1));
        }
        bitset_array[n-1-i] = word;
    }

    return std::make_pair(std::make_pair(bitset_array, n), a.size());


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



#endif //CPU_ENCODERS_AND_DECODERS_H


