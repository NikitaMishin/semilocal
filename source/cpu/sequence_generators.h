//
// Created by nikita on 19.08.2020.
//

#ifndef CPU_SEQUENCE_GENERATORS_H
#define CPU_SEQUENCE_GENERATORS_H


#include <vector>
#include <cstdlib>

/**
 * Generate sequence of int numbers
 * @param size int
 * @param alphabet_size int, default = 4
 * @return vector<int> of size @size
 */
std::vector<int> gen_vector_seq(int size, int alphabet_size = 4) {
    auto v = std::vector<int>();
    v.reserve(size);
    for (int i = 0; i < size; ++i) {
        v.push_back(rand() % alphabet_size);
    }
    return v;
}


#endif //CPU_SEQUENCE_GENERATORS_H