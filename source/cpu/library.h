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


#endif //CPU_LIBRARY_H
