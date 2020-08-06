#ifndef CPU_LIBRARY_H
#define CPU_LIBRARY_H

#include <vector>


#include <iostream>
#include <bitset>
#include <cstring>
#include <map>
#include <unordered_map>
#include "utils.h"
#include "prefix_lcs.h"
#include "transposition_network.h"
#include  <cstdlib>

template<class Input>
int *sticky_braid_sequential(std::vector<Input> &a, std::vector<Input> &b) {
    auto m = a.size();
    auto n = b.size();

    auto size = m + n;
    auto strand_map = new int[size];
    auto reduced_sticky_braid = new int[size];
    for (int i = 0; i < size; ++i) {
        strand_map[i] = i;
    }

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            auto left_edge = m - 1 - i;
            auto top_edge = m + j;
            auto left_strand = strand_map[left_edge];
            auto right_strand = strand_map[top_edge];

            if (a[i] == b[j] || (left_strand > right_strand)) std::swap(strand_map[top_edge], strand_map[left_edge]);

            if (j == n - 1) {
                reduced_sticky_braid[strand_map[left_edge]] = left_edge + n;
            }
            if (i == m - 1) {
                reduced_sticky_braid[strand_map[top_edge]] = top_edge - m;
            }

        }
    }
    return reduced_sticky_braid;
}

/**
 * Assume a <= b
 * @tparam Input
 * @param a
 * @param b
 * @return
 */
template<class Input>
int *sticky_braid_mpi(std::vector<Input> const &a, std::vector<Input> const &b, int threads_num = 1) {

    auto m = a.size();
    auto n = b.size();

    auto size = m + n;
//    int strand_map[100000/2 + 500000*2];// = new int[size];
//    int  reduced_sticky_braid[100000/2 + 500000*2];// = new int[size];
    auto reduced_sticky_braid = new int[size];
    auto strand_map = new int[size];

    auto num_diag = a.size() + b.size() - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);


#pragma omp parallel num_threads(threads_num)  default(none) shared(a, b, strand_map, total_same_length_diag, size, m, n, reduced_sticky_braid)
    {
        int left_edge, top_edge;
        //    init phase
#pragma omp for simd schedule(static) safelen(1)
        for (int k = 0; k < m; ++k) {
            strand_map[k] = k;
        }

#pragma omp for simd schedule(static) safelen(1)
        for (int l = 0; l < n; ++l) {
            strand_map[l + m] = l + m;
        }

        //    phase one
        for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
            left_edge = m - 1 - cur_diag_len;
            top_edge = m;
#pragma omp for simd schedule(static) safelen(1)
            for (int j = 0; j < cur_diag_len + 1; ++j) {
                int left_strand = strand_map[left_edge + j];
                int right_strand = strand_map[top_edge + j];
//                bool r = a[left_edge + j] == b[j] || (!left_strand && right_strand);
                bool r = a[cur_diag_len - j] == b[j] || (left_strand > right_strand);
                if (r) std::swap(strand_map[top_edge + j], strand_map[left_edge + j]);
            }
        }


        for (int j = 0; j < total_same_length_diag; ++j) {
            left_edge = 0;
            top_edge = m + j;
            auto i = m - 1;
#pragma omp for simd schedule(static) safelen(1)
            for (int k = 0; k < m; ++k) {
                auto left_strand = strand_map[left_edge + k];
                auto right_strand = strand_map[top_edge + k];
                bool r = a[i - k] == b[left_edge + j + k] || (left_strand > right_strand);
//                auto r = a[left_edge + k] == b[left_edge + j + k] || (left_strand > right_strand);
                if (r) std::swap(strand_map[top_edge + k], strand_map[left_edge + k]);
            }
        }

////    phase 3
        auto start_j = total_same_length_diag;
        for (int diag_len = m - 2; diag_len >= 0; --diag_len, start_j++) {
            left_edge = 0;
            top_edge = start_j + m;
            auto i = m - 1;
            auto j = start_j;
#pragma omp for simd schedule(static) safelen(1)
            for (int k = 0; k < diag_len + 1; ++k) {
                auto right_strand = strand_map[top_edge + k];
                auto left_strand = strand_map[left_edge + k];
//                auto r = a[left_edge+k] == b[j + k] || (!left_strand && right_strand);
                bool r = a[i - k] == b[j + k] || (left_strand > right_strand);
                if (r) std::swap(strand_map[top_edge + k], strand_map[left_edge + k]);
            }
        }

#pragma omp for simd schedule(static) safelen(1)
        for (int l = 0; l < m; ++l) {
//            reduced_sticky_braid[strand_map[l]] = n + l;
            reduced_sticky_braid[n + l] = strand_map[l]; // seems faster
        }
#pragma omp for simd schedule(static) safelen(1)
        for (int r = m; r < n + m; ++r) {
//            reduced_sticky_braid[strand_map[r]] = r - m;
            reduced_sticky_braid[r - m] = strand_map[r];//seems faster
        }

    }

    delete[] strand_map;
    return reduced_sticky_braid;
}


/**
 *
 * @tparam Input
 * @tparam Output
 * @param a
 * @return
 */
template<class Input,class Output>
std::pair<std::pair<Output *,std::pair<int,int>>, std::unordered_map<Output, Input>> encode(std::vector<Input> const &a) {

    auto mapper_forward = new std::unordered_map<char, Output>();
    auto mapper_reverse = new std::unordered_map<Output, Input>();
    for (int i = 0, encoder = 0; i < a.size(); ++i) {
        if (mapper_forward->count(a[i]) == 0) {
            (*mapper_forward)[a[i]] = Output(encoder);
            (*mapper_reverse)[encoder] = a[i];
            encoder++;
        }
    }
    auto alphabet_size = mapper_reverse->size();
    auto bits_per_symbol = int(std::ceil(log2(alphabet_size)));
    auto shift = bits_per_symbol;
    auto word_size_in_bits = sizeof(Output) * 8;
    auto symbols_in_word = word_size_in_bits / shift;

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

    return std::make_pair(std::make_pair(bitset_array,std::make_pair(n,a.size())), *mapper_reverse);


}




#endif //CPU_LIBRARY_H
