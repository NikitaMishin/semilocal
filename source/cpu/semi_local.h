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
#include <chrono>


void steady_ant_wrapper(AbstractPermutation &p, AbstractPermutation &q, AbstractPermutation &product,
                        distance_unit_monge_product::steady_ant::PrecalcMap &map, int nested_lvls = 3) {


    int nearest_2_degree = pow(2, int(ceil(log2(2 * p.row_size))));
    int total = int(log2(nearest_2_degree)) * nearest_2_degree;
    auto memory_block = new int[p.row_size * 8 + total];
    auto m_new = PermutationPreAllocated(p.row_size, p.col_size, memory_block, memory_block + p.row_size);
    auto n_new = PermutationPreAllocated(p.row_size, p.col_size, memory_block + 2 * p.row_size,
                                         memory_block + 3 * p.row_size);

    copy(p, m_new);
    copy(q, n_new);


    distance_unit_monge_product::steady_ant::steady_ant_with_precalc_and_memory(&m_new, &n_new, memory_block,
                                                                                memory_block + 4 * p.row_size,
                                                                                memory_block + 8 * p.row_size,
                                                                                map, total, nested_lvls);
    copy(m_new, product);

    delete[] memory_block;

}


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
         * Allows get P_{b,a} when you have P_{a,b}
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
                staggered_sticky_multiplication(subtree_l, subtree_r, n, map, product, 0);

                return product;


            }

        }


    }

    /**
     *
     */
    namespace hybrid_approach {
        using namespace distance_unit_monge_product::steady_ant;

        void hybrid(const int *a, int m, const int *b, int n, steady_ant_approach::PrecalcMap &map,
                    AbstractPermutation &perm, int size_bound,
                    int thds_per_combing_algo, int nested_parallel = 0) {


            if (n + m < size_bound) {
                strand_combing_approach::sticky_braid_mpi(perm, a, m, b, n, thds_per_combing_algo);
                return;
            }


            if (n > m) {
                auto n1 = n / 2;
                auto b1 = b;
                auto b2 = b + n1;

                auto subtree_l = Permutation(n1 + m, n1 + m);
                auto subtree_r = Permutation(n - n1 + m, n - n1 + m);

                hybrid(b1, n1, a, m, map, subtree_l, size_bound, thds_per_combing_algo, nested_parallel);
                hybrid(b2, n - n1, a, m, map, subtree_r, size_bound, thds_per_combing_algo, nested_parallel);

                auto product = Permutation(perm.row_size, perm.col_size);

                staggered_sticky_multiplication(&subtree_l, &subtree_r, m, map, &product, 0);
                steady_ant_approach::fill_permutation_ba(&product, &perm, m, n);
            } else {

                auto m1 = m / 2;
                auto a1 = a;
                auto a2 = a + m1;

                auto subtree_l = Permutation(m1 + n, m1 + n);
                auto subtree_r = Permutation(m - m1 + n, m - m1 + n);


                if (nested_parallel > 0) {
#pragma omp parallel num_threads(2)
                    {
#pragma omp single nowait
                        {
#pragma omp task
                            hybrid(a1, m1, b, n, map, subtree_l, size_bound, thds_per_combing_algo,
                                   nested_parallel - 1);
#pragma omp task
                            hybrid(a2, m - m1, b, n, map, subtree_r, size_bound, thds_per_combing_algo,
                                   nested_parallel - 1);
                        }

                    }
#pragma omp taskwait
                } else {
                    hybrid(a1, m1, b, n, map, subtree_l, size_bound, thds_per_combing_algo, nested_parallel);
                    hybrid(a2, m - m1, b, n, map, subtree_r, size_bound, thds_per_combing_algo, nested_parallel);
                }

                staggered_sticky_multiplication(&subtree_l, &subtree_r, n, map, &perm, 0);
            }


        }


        /**
         * Todo not working, check theory
         * @param a
         * @param a_size
         * @param b
         * @param b_size
         * @param permutation
         * @param map
         * @param threads_num
         */
        void first_and_third_phase_merged(const int *a, int a_size, const int *b, int b_size,
                                          AbstractPermutation &permutation,
                                          steady_ant_approach::PrecalcMap &map,
                                          int nested_parall_regions,int threads_num = 1) {
            if (a_size > b_size) {
                return first_and_third_phase_merged(b, b_size, a, a_size, permutation, map, nested_parall_regions,threads_num);
            }

            //assume |a|<=|b|

            auto m = a_size;
            auto n = b_size;


            auto size = m + n;
            int *strand_map = new int[size + 2 * (m - 1)];
            auto third_phase_map_size = m * 2 - 2;
            auto third_phase_map = strand_map + size;

            auto p = Permutation(m + n, m + n);
            auto q = Permutation(third_phase_map_size, third_phase_map_size);

            auto offset = n - (m - 1);

#pragma omp parallel num_threads(threads_num)  default(none) shared(a, b, strand_map, size, m, n, permutation, p, q, offset, third_phase_map, third_phase_map_size)
            {
                int in_third_phase = m - 1;

                //    init phase
#pragma omp for simd schedule(static) nowait
                for (int k = 0; k < (m + n); ++k) {
                    strand_map[k] = k;
                }

                #pragma omp for simd schedule(static) nowait
                for (int k = 0; k < third_phase_map_size; ++k) {
                    if (k < m - 1) {
                        third_phase_map[k] = 2 * k;
                    } else {
                        third_phase_map[k] = (k - (m - 1)) * 2 + 1;
                    }
                }
#pragma omp barrier

                // parallel first and third phase
                for (int diag_number = 0; diag_number < m - 1; ++diag_number) {
                    #pragma omp for simd schedule(static) nowait
                    for (int pos_in_diag = 0; pos_in_diag < m; ++pos_in_diag) {

                        if (pos_in_diag < in_third_phase) {
                            auto top_edge = diag_number + pos_in_diag;
                            auto left_strand = strand_map[pos_in_diag+ m+n];
                            auto top_strand = strand_map[m - 1 + top_edge+m+n];
                            bool r = a[m - 1 - pos_in_diag] == b[offset + top_edge] || (left_strand > top_strand);
                            if (r) std::swap(strand_map[pos_in_diag+m+n], strand_map[m - 1 + top_edge+m+n]);

                        } else {
                            auto top_edge = diag_number + pos_in_diag + 1 - m;
                            auto left_strand = strand_map[pos_in_diag];
                            auto top_strand = strand_map[m + top_edge];
                            bool r = a[m - 1 - pos_in_diag] == b[top_edge] || (left_strand > top_strand);
                            if (r) std::swap(strand_map[pos_in_diag], strand_map[m + top_edge]);
                        }
                    }

                    in_third_phase--;
                    #pragma omp barrier
                }



                // second phase
                for (int j = 0; j < offset; ++j) {
                    auto top_edge = m + j;
                    auto i = m - 1;
#pragma omp for simd schedule(static)
                    for (int k = 0; k < m; ++k) {
                        auto left_strand = strand_map[k];
                        auto right_strand = strand_map[top_edge + k];
                        bool r = a[i - k] == b[j + k] || (left_strand > right_strand);
                        if (r) std::swap(strand_map[top_edge + k], strand_map[k]);
                    }
                }


#pragma omp for simd schedule(static) nowait
                for (int l = 0; l < m; l++) {
                    if (l == m - 1) {
                        p.set_point(strand_map[l], n + l);
                    } else {
                        p.set_point(strand_map[l], l * 2 + offset);
                    }
                }


#pragma omp for simd schedule(static) nowait
                for (int r = m; r < m + n; r++) {
                    if ((r - m) < offset) {
                        p.set_point(strand_map[r], r - m);
                    } else {
                        p.set_point(strand_map[r], (r - m - offset + 1) * 2 + offset - 1);
                    }
                }

#pragma omp for simd schedule(static) nowait
                for (int l = 0; l < m - 1; l++) {
                    q.set_point(third_phase_map[l], m - 1 + l);
                }


#pragma omp for simd schedule(static) nowait
                for (int r = m - 1; r < m + m - 2; r++) {
                    q.set_point(third_phase_map[r], r - (m - 1));
                }

#pragma omp barrier
            }

//            auto begin3 = std::chrono::high_resolution_clock::now(); // or use steady_clock if high_resolution_clock::is_steady is false
            glueing_part_to_whole(&p, &q, map, offset, 1, &permutation,nested_parall_regions);
//            auto time3 = std::chrono::high_resolution_clock::now() - begin3;
//            std::cout << "glue " << std::chrono::duration<double, std::milli>(time3).count()
//                      << std::endl;


            delete[] strand_map;


        }


    }
}

#endif //CPU_SEMI_LOCAL_H
