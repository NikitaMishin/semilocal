//
// Created by nikita on 30.07.2020.
//

#ifndef CPU_TRANSPOSITION_NETWORK_H
#define CPU_TRANSPOSITION_NETWORK_H

#include <vector>

template<class Input>
int prefix_lcs_via_braid_sequential(std::vector<Input> a, std::vector<Input> b) {

    auto m = a.size();
    auto n = b.size();

    auto dis_braid = 0;
    auto size = m + n;
    auto strand_map = new bool[size];

    for (int k = 0; k < m; ++k) {
        strand_map[k] = true;
    }

    for (int l = m; l < m + n; ++l) {
        strand_map[l] = false;
    }


    for (int i = 0; i < m; ++i) {
        for (int j = 0, top_edge = m; j < n; ++j, top_edge++) {
            bool left_strand = strand_map[i];
            bool right_strand = strand_map[top_edge];
            bool r = a[i] == b[j] || (!left_strand && right_strand);
            if (r) std::swap(strand_map[top_edge], strand_map[i]);
//            strand_map[left_edge] = !(!right_strand && (d == left_strand || ((d && !left_strand))));
//            strand_map[left_edge] = !(d || right_strand);
//            strand_map[top_edge] = (d && left_strand) || (left_strand && right_strand);
//            strand_map[top_edge] = (!left_strand && right_strand)? left_strand :    left_strand || right_strand &&d ;
//            strand_map[left_edge] =  right_strand || left_strand || !d ;
        }
    }


    for (int i1 = 0; i1 < m; ++i1) {
        dis_braid += strand_map[i1];
    }


    return m - dis_braid;
}

/**
 * Assume a <= b
 * @tparam Input
 * @param a
 * @param b
 * @return
 */
template<class Input>
int prefix_lcs_via_braid_mpi(std::vector<Input> a, std::vector<Input> b, int threads_num = 1) {

    auto m = a.size();
    auto n = b.size();

    auto dis_braid = 0;
    auto size = m + n;
//    short
    auto strand_map = new int[size];

    auto num_diag = a.size() + b.size() - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);


#pragma omp parallel num_threads(threads_num)  default(none) shared(a, b, strand_map, total_same_length_diag, size, m, n, dis_braid)
    {
        int left_edge, top_edge;
        //    init phase
//#pragma omp for schedule(static)
#pragma omp  for simd schedule(static)
        for (int k = 0; k < m; ++k) {
            strand_map[k] = true;
        }
//        #pragma omp for schedule(static)
#pragma omp for simd schedule(static)
        for (int l = 0; l < n; ++l) {
            strand_map[l + m] = false;
        }

        //    phase one
        for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
            left_edge = m - 1 - cur_diag_len;
            top_edge = m;
//#pragma omp for schedule(static)
#pragma omp for simd schedule(static)
            for (int j = 0; j < cur_diag_len + 1; ++j) {
                auto left_strand = strand_map[left_edge + j];
                auto right_strand = strand_map[top_edge + j];
                auto r = a[cur_diag_len - j] == b[j] || (!left_strand && right_strand);
//                bool r = a[left_edge + j] == b[j] || (!left_strand && right_strand);
                if (r) std::swap(strand_map[top_edge + j], strand_map[left_edge + j]);
            }
        }


        for (int j = 0; j < total_same_length_diag; ++j) {
            left_edge = 0;
            top_edge = m + j;
            auto i = m - 1;
#pragma omp for schedule(static)
//#pragma omp for simd schedule(static)
            for (int k = 0; k < m; ++k) {
                auto left_strand = strand_map[left_edge + k];
                auto right_strand = strand_map[top_edge + k];
//                auto r = a[left_edge + k] == b[left_edge + j + k] || (!left_strand && right_strand);
                bool r = a[i - k] == b[left_edge + j + k] || (!left_strand && right_strand);
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
#pragma omp for schedule(static)
//#pragma omp for simd schedule(static)
            for (int k = 0; k < diag_len + 1; ++k) {
                auto left_strand = strand_map[left_edge + k];
                auto right_strand = strand_map[top_edge + k];
//                auto r = a[left_edge+k] == b[j + k] || (!left_strand && right_strand);
                bool r = a[i - k] == b[j + k] || (!left_strand && right_strand);
                if (r) std::swap(strand_map[top_edge + k], strand_map[left_edge + k]);
            }
        }

#pragma omp for reduction(+:dis_braid) schedule(static)
//#pragma omp for simd reduction(+:dis_braid) schedule(static)
        for (int i1 = 0; i1 < m; ++i1) {
            dis_braid += strand_map[i1];
        }

    }

    delete[] strand_map;


    return m - dis_braid;
}


#endif //CPU_TRANSPOSITION_NETWORK_H
