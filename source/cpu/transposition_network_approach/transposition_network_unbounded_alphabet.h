//
// Created by nikita on 19.08.2020.
//

#ifndef CPU_TRANSPOSITION_NETWORK_UNBOUNDED_ALPHABET_H
#define CPU_TRANSPOSITION_NETWORK_UNBOUNDED_ALPHABET_H

#include <vector>
#include <cmath>

/**
 * @tparam Input
 * @tparam StrandHolder
 * @param a
 * @param b
 * @return
 */
template<class Input,class StrandHolder>
int prefix_lcs_via_braid_mpi(std::vector<Input> const &a, std::vector<Input> const &b, int threads_num = 1) {

    if(a.size()>b.size()) return prefix_lcs_via_braid_mpi(b,a,threads_num);

    auto m = a.size();
    auto n = b.size();

    auto dis_braid = 0;
    auto size = m + n;

    auto strand_map = new StrandHolder[size];
    auto num_diag = a.size() + b.size() - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);

    #pragma omp parallel num_threads(threads_num)  default(none) shared(a, b, strand_map, total_same_length_diag, size, m, n, dis_braid)
    {
        int left_edge, top_edge;
        //    init phase
        #pragma omp  for simd schedule(static)
        for (int k = 0; k < m; ++k) {
            strand_map[k] = StrandHolder(true);
        }

        #pragma omp for simd schedule(static)
        for (int l = 0; l < n; ++l) {
            strand_map[l + m] = StrandHolder(false);
        }

        //    phase one
        for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
            left_edge = m - 1 - cur_diag_len;
            top_edge = m;
            #pragma omp for simd schedule(static)
            for (int j = 0; j < cur_diag_len + 1; ++j) {
                StrandHolder left_strand = strand_map[left_edge + j];
                StrandHolder right_strand = strand_map[top_edge + j];
                if ((a[cur_diag_len - j] == b[j]) || (left_strand > right_strand)) std::swap(strand_map[top_edge + j], strand_map[left_edge + j]);
            }
        }


        for (int j = 0; j < total_same_length_diag; ++j) {
            top_edge = m + j;
            auto i = m - 1;
            #pragma omp for simd schedule(static)
            for (int k = 0; k < m; ++k) {
                StrandHolder left_strand = strand_map[k];
                StrandHolder right_strand = strand_map[top_edge + k];
                if ((a[i - k] == b[j + k]) || (left_strand > right_strand)) std::swap(strand_map[top_edge + k], strand_map[k]);
            }
        }

////    phase 3
        auto start_j = total_same_length_diag;
        for (int diag_len = m - 2; diag_len >= 0; --diag_len, start_j++) {
            top_edge = start_j + m;
            auto i = m - 1;
            auto j = start_j;
            #pragma omp for simd schedule(static)
            for (int k = 0; k < diag_len + 1; ++k) {
                StrandHolder left_strand = strand_map[k];
                StrandHolder right_strand = strand_map[top_edge + k];
                if ((a[i - k] == b[j + k]) || (left_strand > right_strand)) std::swap(strand_map[top_edge + k], strand_map[k]);
            }
        }

        #pragma omp for simd reduction(+:dis_braid) schedule(static)
        for (int i1 = 0; i1 < m; ++i1) {
            dis_braid += strand_map[i1];
        }
    }
    delete[] strand_map;

    return m - dis_braid;
}



template<class Input,class StrandHolder>
int prefix_lcs_via_braid_sequential(std::vector<Input> const &a, std::vector<Input> const &b) {

    auto m = a.size();
    auto n = b.size();

    auto dis_braid = 0;
    auto size = m + n;
    auto strand_map = new StrandHolder[size];

    for (int k = 0; k < m; ++k) {
        strand_map[k] = StrandHolder(true);
    }

    for (int l = m; l < m + n; ++l) {
        strand_map[l] = StrandHolder(false);
    }


    for (int i = 0; i < m; ++i) {
        for (int j = 0, top_edge = m; j < n; ++j, top_edge++) {
            StrandHolder left_strand = strand_map[i];
            StrandHolder right_strand = strand_map[top_edge];
            if ((a[i] == b[j]) || (!left_strand && right_strand)) std::swap(strand_map[top_edge], strand_map[i]);
        }
    }

    for (int i1 = 0; i1 < m; ++i1) {
        dis_braid += strand_map[i1];
    }


    return m - dis_braid;
}









#endif //CPU_TRANSPOSITION_NETWORK_UNBOUNDED_ALPHABET_H




