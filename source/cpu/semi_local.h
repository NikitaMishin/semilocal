#ifndef CPU_SEMI_LOCAL_H
#define CPU_SEMI_LOCAL_H

#include <vector>


#include <iostream>
#include <bitset>
#include <cstring>
#include <map>
#include <unordered_map>
#include "unit_monge_mult/steady_ant.h"
#include  <cstdlib>
#include <chrono>
#include "predefined_types.h"
#include "deque"

struct{
    int fst_1;
    int snd_1;
    int fst_2;
    int snd_2;
    bool bit;

} tripple;


void steady_ant_wrapper(AbstractPermutation &p, AbstractPermutation &q, AbstractPermutation &product,
                        PrecalcMap &map, int nested_lvls = 3) {


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

    namespace _details {
        using namespace distance_unit_monge_product::steady_ant;

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
        anti_diagonal_computation(Input *strand_map, const Input *a, const Input *b, int upper_bound, int left_edge,
                                  int top_edge, int offset_a, int offset_b) {

            #pragma omp for simd schedule(static) aligned(a, b, strand_map:sizeof(int)*8) nowait
            for (int k = 0; k < upper_bound; ++k) {

                auto left_strand = strand_map[left_edge + k];
                auto right_strand = strand_map[top_edge + k];

                bool r = a[offset_a + k] == b[offset_b + k] || (left_strand > right_strand);

                if (WithIf) {
                    if (r) {
                        strand_map[top_edge + k] = left_strand;
                        strand_map[left_edge + k] = right_strand;
                    }
                } else {
                    strand_map[left_edge + k] = (left_strand & (r - 1)) | ((-r) & right_strand);
                    strand_map[top_edge + k] = (right_strand & (r - 1)) | ((-r) & left_strand);
                }
            }

            if (withWait){
                #pragma omp barrier
            }
        }




        template<class Input>
        inline void initialization(Input *strand_map, int m, int n) {
#pragma omp for simd schedule(static)
            for (int k = 0; k < m; ++k) {
                strand_map[k] = k;
            }

#pragma omp for simd schedule(static)
            for (int l = 0; l < n; ++l) {
                strand_map[l + m] = l + m;
            }

        }

        template<class Input>
        inline void
        construct_permutation(AbstractPermutation &matrix, Input *strand_map, bool is_reverse, int m, int n) {
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

        template<class Input>
        inline void fill_a_reverse(const Input *a, Input *a_reverse, int m) {
        #pragma omp  for simd schedule(static)
            for (int i = 0; i < m; ++i) {
                a_reverse[i] = a[m - 1 - i];
            }
        }


    }


    template<class Input, bool WithIf>
    void
    sticky_braid_sequential(AbstractPermutation &permutation, const Input *a, int a_size, const Input *b, int b_size) {

        auto m = a_size;
        auto n = b_size;

        auto top_strands = new Input[n];
        auto left_strands = new Input[m];

        // init phase
        for (int i = 0; i < m; ++i) {
            left_strands[i] = i;
        }
        for (int i = 0; i < n; ++i) {
            top_strands[i] = i + m;
        }

        for (int i = 0; i < m; ++i) {
            auto left_edge = m - 1 - i;
            auto left_strand = left_strands[left_edge];
            auto a_symbol = a[i];
            int right_strand;
            for (int j = 0; j < n - 1; ++j) {
                right_strand = top_strands[j];
                auto r = a_symbol == b[j] || (left_strand > right_strand);

                if (WithIf) {
                    if (r) {
                        top_strands[j] = left_strand;
                        left_strand = right_strand;
                    }
                } else {

                    top_strands[j] = (right_strand & (r - 1)) | ((-r) & left_strand);
                    left_strand = (left_strand & (r - 1)) | ((-r) & right_strand);

                }

            }

            right_strand = top_strands[n - 1];
            auto r = a_symbol == b[n - 1] || (left_strand > right_strand);
            left_strands[left_edge] = (left_strand & (r - 1)) | ((-r) & right_strand);

            if (WithIf) {
                if (r) top_strands[n - 1] = left_strand;
            } else {
                top_strands[n - 1] = (right_strand & (r - 1)) | ((-r) & left_strand);
            }

        }


        // permutation construction phase
        for (int l = 0; l < m; l++) permutation.set_point(left_strands[l], n + l);

        for (int r = m; r < m + n; r++) permutation.set_point(top_strands[r - m], r - m);

        delete[] left_strands;
        delete[] top_strands;
    }


    /**
     *
     * @tparam Input
     * @tparam WithIf
     * @param matrix
     * @param a
     * @param a_size
     * @param b
     * @param b_size
     * @param threads_num
     * @param is_reverse
     */
    template<class Input, bool WithIf,bool WithWait>
    void sticky_braid_mpi(AbstractPermutation &matrix, const Input *a, int a_size, const Input *b, int b_size,
                          int threads_num = 1, bool is_reverse = false) {
        using namespace _details;

        if (a_size > b_size) {
            sticky_braid_mpi<Input, WithIf,WithWait>(matrix, b, b_size, a, a_size, threads_num, !is_reverse);
            return;
        }

        auto m = a_size;
        auto n = b_size;


        auto size = m + n;
        Input *strand_map = new Input[size];

        auto num_diag = m + n - 1;
        auto total_same_length_diag = num_diag - (m - 1) - (m - 1);
        Input *a_reverse = new Input[m];

#pragma omp parallel num_threads(threads_num)  default(none) shared(a_reverse, a, b, is_reverse, strand_map, matrix, total_same_length_diag, size, m, n)
        {


            int left_edge, top_edge;
            //    init phase
            initialization(strand_map, m, n);
            fill_a_reverse(a, a_reverse, m);

            //    phase one
            top_edge = m;
            left_edge = m - 1;
            for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
                anti_diagonal_computation<Input, WithIf,WithWait>(strand_map, a_reverse, b, cur_diag_len + 1, left_edge,
                                                         top_edge, left_edge, 0);
                left_edge--;
            }

            //phase 2
            top_edge = m;
            for (int j = 0; j < total_same_length_diag; ++j) {
                anti_diagonal_computation<Input, WithIf,WithWait>(strand_map, a_reverse, b, m, 0, top_edge, 0, j);
                top_edge++;
            }

            //// phase 3
            auto start_j = total_same_length_diag;
            top_edge = start_j + m;

            for (int diag_len = m - 2; diag_len >= 0; --diag_len, start_j++) {
                anti_diagonal_computation<Input, WithIf,WithWait>(strand_map, a_reverse, b, diag_len + 1, 0, top_edge, 0,
                                                         start_j);
                top_edge++;
            }

            construct_permutation(matrix, strand_map, is_reverse, m, n);
        }

        delete[] a_reverse;
        delete[] strand_map;
    }

    template<class Input, bool WithIf>
    void
    first_and_third_phase_combined(AbstractPermutation &matrix, const Input *a, int a_size, const Input *b,
                                   int b_size, PrecalcMap &map, int nested_parall_regions = 0, int threads_num = 1) {
        using namespace _details;

        if (a_size > b_size) {
            auto p = Permutation(a_size + b_size, a_size + b_size);
            first_and_third_phase_combined<Input, WithIf>(p, b, b_size, a, a_size, map, nested_parall_regions, threads_num);
            fill_permutation_ba(&p, &matrix, a_size, b_size);
            return;
        }

        //assume |a|<=|b|

        auto m = a_size;
        auto n = b_size;

        auto size = m + n;
        Input *strand_map = new Input[size + 2 * (m - 1)];
        auto third_phase_map_size = m * 2 - 2;
        auto third_phase_map = strand_map + size;

        auto p = Permutation(m + n, m + n);
        auto q = Permutation(third_phase_map_size, third_phase_map_size);

        auto offset = n - (m - 1);
        Input *a_reverse = new Input[m];

#pragma omp parallel num_threads(threads_num)  default(none) shared(a, a_reverse, b, strand_map, size, m, n, matrix, p, q, offset, third_phase_map, third_phase_map_size)
        {

            fill_a_reverse(a, a_reverse, m);
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


            for (int diag_number = 0; diag_number < m - 1; ++diag_number) {


                #pragma omp for simd schedule(static) nowait
                for (int pos_in_diag = 0; pos_in_diag < in_third_phase; ++pos_in_diag) {

                    auto top_edge = diag_number + pos_in_diag;
                    auto left_strand = third_phase_map[pos_in_diag];
                    auto top_strand = third_phase_map[m - 1 + top_edge];
                    bool r = a_reverse[pos_in_diag] == b[offset + top_edge] || (left_strand > top_strand);

                    if (WithIf) {
                        if (r) std::swap(third_phase_map[pos_in_diag], third_phase_map[m - 1 + top_edge]);
                    } else {
                        third_phase_map[pos_in_diag] = (left_strand & (r - 1))  | ((-r) &  top_strand);
                        third_phase_map[m - 1 + top_edge]  = (top_strand & (r - 1)) | ((-r) & left_strand );
                    }


                }

                #pragma omp for simd schedule(static)
                for (int pos_in_diag = in_third_phase; pos_in_diag < m; ++pos_in_diag) {
                    auto top_edge = diag_number + pos_in_diag + 1 - m;
                    auto left_strand = strand_map[pos_in_diag];
                    auto top_strand = strand_map[m + top_edge];
                    bool r = a_reverse[pos_in_diag] == b[top_edge] || (left_strand > top_strand);

                    if (WithIf) {
                        if (r) if (r) std::swap(strand_map[pos_in_diag], strand_map[m + top_edge]);
                    } else {
                        strand_map[pos_in_diag] = (left_strand & (r - 1))  | ((-r) &  top_strand);
                        strand_map[m + top_edge]  = (top_strand & (r - 1)) | ((-r) & left_strand );
                    }


                }
                in_third_phase--;
            }

            //phase 2
            auto top_edge = m;
            for (int j = 0; j < offset; ++j) {
                anti_diagonal_computation<Input, WithIf>(strand_map, a_reverse, b, m, 0, top_edge, 0, j);
                top_edge++;
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

        glueing_part_to_whole(&p, &q, map, offset, 1, &matrix, nested_parall_regions);

        delete[] strand_map;
        delete[] a_reverse;
    }





    template<class Input, bool WithIf, bool WithWait, bool UseSumBound>
    void hybrid(
            AbstractPermutation &perm, const Input *a, int m, const Input *b, int n,
            PrecalcMap &map, int thds_per_combing_algo,
            int braid_mul_parall_depth, int depth, int sum_bound, int parallel_depth) {

        using namespace _details;


        if(UseSumBound) {
            if (m + n <= sum_bound) {
                sticky_braid_mpi<Input, WithIf,WithWait>(perm, a, m, b, n, thds_per_combing_algo);
                return;
            }

        } else {
            //base case
            if (depth <= 0) {
                sticky_braid_mpi<Input, WithIf,WithWait>(perm, a, m, b, n, thds_per_combing_algo);
                return;
            }
        }


        if (n > m) {
            auto n1 = n / 2;
            auto b1 = b;
            auto b2 = b + n1;

            auto subtree_l = Permutation(n1 + m, n1 + m);
            auto subtree_r = Permutation(n - n1 + m, n - n1 + m);

            if (parallel_depth > 0) {
#pragma omp parallel num_threads(2)
                {
#pragma omp single nowait
                    {
#pragma omp task
                        hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_l, b1, n1, a, m, map, thds_per_combing_algo,
                                              braid_mul_parall_depth, depth - 1, sum_bound, parallel_depth - 1);

#pragma omp task
                        hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_r, b2, n - n1, a, m, map, thds_per_combing_algo,
                                              braid_mul_parall_depth, depth - 1, sum_bound, parallel_depth - 1);
                    }

                }
#pragma omp taskwait
            } else {
                hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_l, b1, n1, a, m, map, thds_per_combing_algo, braid_mul_parall_depth,
                                      depth - 1, sum_bound, parallel_depth - 1);
                hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_r, b2, n - n1, a, m, map, thds_per_combing_algo,
                                      braid_mul_parall_depth, depth - 1, sum_bound, parallel_depth - 1);
            }

            auto product = Permutation(perm.row_size, perm.col_size);

            staggered_sticky_multiplication(&subtree_l, &subtree_r, m, map, &product, braid_mul_parall_depth);
            fill_permutation_ba(&product, &perm, m, n);
        } else {

            auto m1 = m / 2;
            auto a1 = a;
            auto a2 = a + m1;

            auto subtree_l = Permutation(m1 + n, m1 + n);
            auto subtree_r = Permutation(m - m1 + n, m - m1 + n);


            if (parallel_depth > 0) {
#pragma omp parallel num_threads(2)
                {
#pragma omp single nowait
                    {
#pragma omp task
                        hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_l, a1, m1, b, n, map, thds_per_combing_algo,
                                              braid_mul_parall_depth, depth - 1, sum_bound, parallel_depth - 1);
#pragma omp task
                        hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_r, a2, m - m1, b, n, map, thds_per_combing_algo,
                                              braid_mul_parall_depth, depth - 1, sum_bound, parallel_depth - 1);
                    }
                }
#pragma omp taskwait
            } else {
                hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_l, a1, m1, b, n, map, thds_per_combing_algo, braid_mul_parall_depth,
                                      depth - 1, sum_bound, parallel_depth - 1);
                hybrid<Input, WithIf,WithWait, UseSumBound>(subtree_r, a2, m - m1, b, n, map, thds_per_combing_algo,
                                      braid_mul_parall_depth, depth - 1, sum_bound, parallel_depth - 1);
            }

            staggered_sticky_multiplication(&subtree_l, &subtree_r, n, map, &perm, braid_mul_parall_depth);
        }


    }


    template<class Input, bool WithIf, bool WithWait, bool UseSumBound>
    void hybrid_non_rec(
            AbstractPermutation &perm, const Input *a, int m, const Input *b, int n,
            PrecalcMap &map  ) {

        using namespace _details;
        // TDOO implement buttom-up appaorch




    }


}

#endif //CPU_SEMI_LOCAL_H
