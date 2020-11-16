#ifndef CPU_SEMI_LOCAL_H
#define CPU_SEMI_LOCAL_H

#include <vector>


#include <iostream>
#include <bitset>
#include <cstring>
#include <map>
#include <unordered_map>
#include "transposition_network_approach/transposition_network_4symbol_alphabet_bit.h"
#include "unit_monge_mult/steady_ant.h"
#include  <cstdlib>


namespace semi_local {


    /**
     *
     */
    namespace strand_combing_approach {

        void
        sticky_braid_sequential(AbstractPermutation &permutation, const int *a, int a_size, const int *b, int b_size) {
            auto m = a_size;
            auto n = b_size;

            auto size = m + n;
            int *strand_map = new int[size];
            for (int i = 0; i < size; ++i) {
                strand_map[i] = i;
            }

            for (int i = 0; i < m; ++i) {
                for (int j = 0; j < n; ++j) {
                    auto left_edge = m - 1 - i;
                    auto top_edge = m + j;
                    auto left_strand = strand_map[left_edge];
                    auto right_strand = strand_map[top_edge];

                    if (a[i] == b[j] || (left_strand > right_strand))
                        std::swap(strand_map[top_edge], strand_map[left_edge]);

                    if (j == n - 1) {
                        permutation.set_point(strand_map[left_edge], left_edge + n);
                    }
                    if (i == m - 1) {
                        permutation.set_point(strand_map[top_edge], top_edge - m);
                    }

                }
            }

            delete[] strand_map;

        }

        void sticky_braid_mpi(AbstractPermutation &matrix, const int *a, int a_size, const int *b, int b_size,
                              int threads_num = 1, bool is_reverse = false) {

            if (a_size > b_size) {
                sticky_braid_mpi(matrix, b, b_size, a, a_size, threads_num, !is_reverse);

                return;
            }

            auto m = a_size;
            auto n = b_size;

            auto size = m + n;
            int *strand_map = new int[size];

            auto num_diag = m + n - 1;
            auto total_same_length_diag = num_diag - (m - 1) - (m - 1);


#pragma omp parallel num_threads(threads_num)  default(none) shared(a, b, is_reverse, strand_map, matrix, total_same_length_diag, size, m, n)
            {
                int left_edge, top_edge;
                //    init phase
#pragma omp for simd schedule(static)
                for (int k = 0; k < m; ++k) {
                    strand_map[k] = k;
                }

#pragma omp for simd schedule(static)
                for (int l = 0; l < n; ++l) {
                    strand_map[l + m] = l + m;
                }

                //    phase one
                for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
                    left_edge = m - 1 - cur_diag_len;
                    top_edge = m;
#pragma omp for simd schedule(static)
                    for (int j = 0; j < cur_diag_len + 1; ++j) {
                        auto left_strand = strand_map[left_edge + j];
                        auto right_strand = strand_map[top_edge + j];
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
                        auto left_strand = strand_map[left_edge + k];
                        auto right_strand = strand_map[top_edge + k];
                        bool r = a[i - k] == b[left_edge + j + k] || (left_strand > right_strand);
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
                        auto right_strand = strand_map[top_edge + k];
                        auto left_strand = strand_map[left_edge + k];

                        bool r = a[i - k] == b[j + k] || (left_strand > right_strand);
                        if (r) std::swap(strand_map[top_edge + k], strand_map[left_edge + k]);
                    }
                }

                if (!is_reverse) {
#pragma omp for simd schedule(static)
                    for (int r = m; r < m + n; r++) {
                        matrix.set_point(strand_map[r], r - m);
                    }
#pragma omp for simd schedule(static)
                    for (int l = 0; l < m; l++) {
                        matrix.set_point(strand_map[l], n + l);
                    }


                } else {

#pragma omp for simd schedule(static)
                    for (int r = m; r < m + n; r++) {
                        matrix.set_point(n + m - 1 - strand_map[r], n + m - 1 - (r - m));
                    }
#pragma omp for simd schedule(static)
                    for (int l = 0; l < m; l++) {
                        matrix.set_point(n + m - 1 - strand_map[l], n + m - 1 - (n + l));
                    }

                }


            }

            delete[] strand_map;
        }


        template<class Input, class StrandHolder>
        StrandHolder *sticky_braid_sequential_without_if(std::vector<Input> &a, std::vector<Input> &b) {
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
                    int c = a[i] == b[j] || (left_strand > right_strand);

                    strand_map[top_edge] = (1 - c) * right_strand + c * left_strand;
                    strand_map[left_edge] = (1 - c) * left_strand + c * right_strand;


                }
            }

            for (int l = 0; l < m; ++l) {
//            reduced_sticky_braid[strand_map[l]] = n + l;
                reduced_sticky_braid[n + l] = strand_map[l]; // seems faster
            }

            for (int r = m; r < n + m; ++r) {
//            reduced_sticky_braid[strand_map[r]] = r - m;
                reduced_sticky_braid[r - m] = strand_map[r];//seems faster
            }


            delete[] strand_map;

            return reduced_sticky_braid;
        }


    }

    /**
     *
     */
    namespace steady_ant_approach {
        using namespace distance_unit_monge_product::steady_ant;

        /**
        * see theorem 5.21
         * Allows get P_{a,b} when you have P_{b,a}
        */
        void fill_permutation_ba(AbstractPermutation *ab, AbstractPermutation *ba, int m, int n) {
            ba->unset_all();
            for (int i = 0; i < ab->row_size; ++i) {
                auto col = ab->get_col_by_row(i);
                if (col != NOPOINT) ba->set_point(n + m - 1 - i, m + n - 1 - col);
            }
        }

        AbstractPermutation *get_semi_local_kernel(int *a, int m, int *b, int n, PrecalcMap &map) {
            if (m == 1 && n == 1) {
                auto p = new Permutation(2, 2);
                if (a[0] == b[0]) {
                    p->set_point(0, 0);
                    p->set_point(1, 1);
                } else {
                    p->set_point(0, 1);
                    p->set_point(1, 0);
                }
                return p;

            }

            if (n > m) {
                auto n1 = n / 2;
                auto b1 = b;
                auto b2 = b + n1;

                auto subtree_l = get_semi_local_kernel(b1, n1, a, m, map);
                auto subtree_r = get_semi_local_kernel(b2, n - n1, a, m, map);
                auto product = new Permutation(subtree_l->row_size + subtree_r->row_size - m,
                                               subtree_l->col_size + subtree_r->col_size - m);

                auto product_t = new Permutation(subtree_l->col_size + subtree_r->col_size - m,
                                                 subtree_l->row_size + subtree_r->row_size - m);


                staggered_sticky_multiplication(subtree_l, subtree_r, m, map, product);
                fill_permutation_ba(product, product_t, m, n);
                return product_t;
            } else {
                auto m1 = m / 2;
                auto a1 = a;
                auto a2 = a + m1;

                auto subtree_l = get_semi_local_kernel(a1, m1, b, n, map);
                auto subtree_r = get_semi_local_kernel(a2, m - m1, b, n, map);
                auto product = new Permutation(subtree_l->row_size + subtree_r->row_size - n,
                                               subtree_l->col_size + subtree_r->col_size - n);
                staggered_sticky_multiplication(subtree_l, subtree_r, n, map, product);

                return product;


            }

        }


    }

    /**
     *
     */
    namespace hybrid_approach {

        




    }


}


#endif //CPU_SEMI_LOCAL_H
