#ifndef CPU_SEMI_LOCAL_H
#define CPU_SEMI_LOCAL_H

#include <vector>


#include <iostream>
#include <bitset>
#include <cstring>
#include <map>
#include <unordered_map>
#include "transposition_network_approach/transposition_network_4symbol_alphabet_bit.h"
#include  <cstdlib>


template<class Input, class StrandHolder>
StrandHolder *sticky_braid_sequential(std::vector<Input> &a, std::vector<Input> &b) {
    auto m = a.size();
    auto n = b.size();

    auto size = m + n;
    StrandHolder *strand_map = new StrandHolder[size];
    StrandHolder *reduced_sticky_braid = new StrandHolder[size];
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
//                reduced_sticky_braid[strand_map[left_edge]] = left_edge + n;
                reduced_sticky_braid[left_edge + n] = strand_map[left_edge];
            }
            if (i == m - 1) {
//                reduced_sticky_braid[strand_map[top_edge]] = top_edge - m;
                reduced_sticky_braid[top_edge - m] = strand_map[top_edge];
            }

        }
    }

    delete [] strand_map;

    return reduced_sticky_braid;
}

/**
 * @tparam Input
 * @param a
 * @param b
 * @return
 */
template<class Input, class StrandHolder>
StrandHolder *sticky_braid_mpi(std::vector<Input> const &a, std::vector<Input> const &b, int threads_num = 1) {

    if (a.size() > b.size()) return sticky_braid_mpi<Input, StrandHolder>(b, a, threads_num);

    auto m = a.size();
    auto n = b.size();

    auto size = m + n;
    StrandHolder *reduced_sticky_braid = new StrandHolder[size];
    StrandHolder *strand_map = new StrandHolder[size];

    auto num_diag = a.size() + b.size() - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);


#pragma omp parallel num_threads(threads_num)  default(none) shared(a, b, strand_map, total_same_length_diag, size, m, n, reduced_sticky_braid)
    {
        int left_edge, top_edge;
        //    init phase
#pragma omp for simd schedule(static)
        for (int k = 0; k < m; ++k) {
            strand_map[k] = StrandHolder(k);
        }

#pragma omp for simd schedule(static)
        for (int l = 0; l < n; ++l) {
            strand_map[l + m] = StrandHolder(l + m);
        }

        //    phase one
        for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
            left_edge = m - 1 - cur_diag_len;
            top_edge = m;
#pragma omp for simd schedule(static)
            for (int j = 0; j < cur_diag_len + 1; ++j) {
                StrandHolder left_strand = strand_map[left_edge + j];
                StrandHolder right_strand = strand_map[top_edge + j];
//                bool r = a[left_edge + j] == b[j] || (!left_strand && right_strand);
                bool r = a[cur_diag_len - j] == b[j] || (left_strand > right_strand);
                if (r) std::swap(strand_map[top_edge + j], strand_map[left_edge + j]);
            }
        }


        for (int j = 0; j < total_same_length_diag; ++j) {
            left_edge = 0;
            top_edge = m + j;
            auto i = m - 1;
#pragma omp for simd schedule(static)
            for (int k = 0; k < m; ++k) {
                StrandHolder left_strand = strand_map[left_edge + k];
                StrandHolder right_strand = strand_map[top_edge + k];
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
#pragma omp for simd schedule(static)
            for (int k = 0; k < diag_len + 1; ++k) {
                StrandHolder right_strand = strand_map[top_edge + k];
                StrandHolder left_strand = strand_map[left_edge + k];
//                auto r = a[left_edge+k] == b[j + k] || (!left_strand && right_strand);
                bool r = a[i - k] == b[j + k] || (left_strand > right_strand);
                if (r) std::swap(strand_map[top_edge + k], strand_map[left_edge + k]);
            }
        }

#pragma omp for simd schedule(static)
        for (int l = 0; l < m; ++l) {
//            reduced_sticky_braid[strand_map[l]] = n + l;
            reduced_sticky_braid[n + l] = strand_map[l]; // seems faster
        }

#pragma omp for simd schedule(static)
        for (int r = m; r < n + m; ++r) {
//            reduced_sticky_braid[strand_map[r]] = r - m;
            reduced_sticky_braid[r - m] = strand_map[r];//seems faster
        }

    }

    delete[] strand_map;
    return reduced_sticky_braid;
}


#endif //CPU_SEMI_LOCAL_H
