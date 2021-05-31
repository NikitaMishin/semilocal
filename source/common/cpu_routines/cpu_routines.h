#ifndef SEMI_CPU_ROUTINES_H
#define SEMI_CPU_ROUTINES_H

#include "monge/matrices.h"




/**
 *
 * @tparam Input
 * @tparam WithIf weather or not to use approach with if rather then  the branchless one @param strand_map
 * @param a
 * @param b
 * @param upper_bound
 * @param left_edge
 * @param top_edge
 * @param offset_a
 * @param offset_b
 */
template<class Input, bool WithIf, bool withWait>
inline void
anti_diagonal_computation_cpu(Input *l_strands, Input *t_strands, const Input *a, const Input *b, int upper_bound,
                              int offset_a, int offset_b) {

    #pragma omp for simd schedule(static) aligned(a, b, strand_map:sizeof(Input)*8) nowait
    for (int k = 0; k < upper_bound; ++k) {

        auto left_strand = l_strands[offset_a + k];
        auto right_strand = t_strands[offset_b + k];

        auto r = (a[offset_a + k] == b[offset_b + k]) || (left_strand > right_strand);

        if (WithIf) {
            if (r) {
                t_strands[offset_b + k] = left_strand;
                l_strands[offset_a + k] = right_strand;
            }
        } else {
            auto r_minus = (r - 1);
            auto minus_r = -r;
            auto l_new = (left_strand & r_minus) | (minus_r & right_strand);
            auto r_new = (right_strand & r_minus) | (minus_r & left_strand);

            l_strands[offset_a + k] = l_new;
            t_strands[offset_b + k] = r_new;
        }
    }

    if (withWait) {
        #pragma omp barrier
    }
}

template<class Input>
inline void initialization_cpu(Input *l_strands, Input * t_strands, int m, int n) {
    #pragma omp for simd schedule(static)
    for (int k = 0; k < m; ++k) {
        l_strands[k] = k;
    }

    #pragma omp for simd schedule(static)
    for (int l = 0; l < n; ++l) {
        t_strands[l] = l + m;
    }
}

template<class Input>
inline void
construct_permutation_cpu(AbstractPermutation &matrix, Input *l_strands, Input *t_strands, bool is_reverse, int m, int n) {
    if (!is_reverse) {
        #pragma omp for simd schedule(static)
        for (int r = 0; r < n; r++) {
            matrix.set_point(t_strands[r], r);
        }

        #pragma omp for simd schedule(static)
        for (int l = 0; l < m; l++) {
            matrix.set_point(l_strands[l], n + l);
        }

    } else {

        #pragma omp for simd schedule(static)
        for (int r = 0; r < n; r++) {
            matrix.set_point(n + m - 1 - t_strands[r], n + m - 1 - r);
        }
        #pragma omp for simd schedule(static)
        for (int l = 0; l < m; l++) {
            matrix.set_point(n + m - 1 - l_strands[l], n + m - 1 - (n + l));
        }

    }

}

template<class Input>
inline void fill_a_reverse_cpu(const Input *a, Input *a_reverse, int m) {
    #pragma omp  for simd schedule(static)
    for (int i = 0; i < m; ++i) {
        a_reverse[i] = a[m - 1 - i];
    }
}


void fill_permutation_ba(AbstractPermutation *ab, AbstractPermutation *ba, int m, int n) {
    ba->unset_all();
    for (int i = 0; i < ab->row_size; ++i) {
        auto col = ab->get_col_by_row(i);
        if (col != NOPOINT) ba->set_point(n + m - 1 - i, m + n - 1 - col);
    }
}

#endif //SEMI_CPU_ROUTINES_H
