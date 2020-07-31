//
// Created by nikita on 30.07.2020.
//

#ifndef CPU_UTILS_H
#define CPU_UTILS_H




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
