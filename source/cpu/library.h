#ifndef CPU_LIBRARY_H
#define CPU_LIBRARY_H

#include <vector>


#include <iostream>
#include <bitset>
#include <cstring>
#include "utils.h"
#include "prefix_lcs.h"
#include "transposition_network.h"


template<class Input>
int *sticky_braid_sequential(std::vector<Input> a, std::vector<Input> b) {
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
int * sticky_braid_sequential_skewed(std::vector<Input> a, std::vector<Input> b) {

//    10
//    010
    auto m = a.size();
    auto n = b.size();

    auto size = m + n;
    auto reduced_sticky_braid = new int[size];
    auto strand_map = new int[size];

    auto num_diag = a.size() + b.size() - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);
    int left_edge, top_edge;

    //    init phase
    for (int k = 0; k < m+n; ++k) {
        strand_map[k] = k;
    }

    //    phase one
    for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
        left_edge = m - 1 - cur_diag_len;
        top_edge = m;
        for (int j = 0, i = cur_diag_len; j <= cur_diag_len; ++j, left_edge++, top_edge++, i--) {
            auto left_strand = strand_map[left_edge];
            auto right_strand = strand_map[top_edge];
            bool r = a[i] == b[j] || (left_strand > right_strand);
            if (r) std::swap(strand_map[top_edge], strand_map[left_edge]);
        }


    }

    // phase two
    // equals
    for (int j = 0; j < total_same_length_diag; ++j) {
        left_edge = 0;
        top_edge = m + j;
        auto i = m - 1;
        for (int k = 0; k < m; ++k, top_edge++, left_edge++, i--) {
            auto left_strand = strand_map[left_edge];
            auto right_strand = strand_map[top_edge];
            bool r = a[i] == b[left_edge + j] || (left_strand > right_strand);
            if (r) std::swap(strand_map[top_edge], strand_map[left_edge]);

        }
    }

////    phase 3
    auto start_j = total_same_length_diag;
    for (int diag_len = m - 2; diag_len >= 0; --diag_len, start_j++) {
        left_edge = 0;
        top_edge = start_j + m;
        auto i = m - 1;
        auto j = start_j;
        for (int k = 0; k <= diag_len; ++k, left_edge++, top_edge++, i--, j++) {
            auto left_strand = strand_map[left_edge];
            auto right_strand = strand_map[top_edge];
            bool r = a[i] == b[j] || (left_strand > right_strand);
            if (r) std::swap(strand_map[top_edge], strand_map[left_edge]);

        }

    }

    for (int l = 0; l < m; ++l) {
        reduced_sticky_braid[strand_map[l]] =  n + l;
//        reduced_sticky_braid[n+l] =  strand_map[l];
    }
    for (int r = m; r < n+m ; ++r) {
        reduced_sticky_braid[strand_map[r]] = r - m;
//        reduced_sticky_braid[r-m] = strand_map[r];
    }


    return reduced_sticky_braid;
}


#endif //CPU_LIBRARY_H
