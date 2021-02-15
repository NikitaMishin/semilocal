//
// Created by nikita on 19.08.2020.
//

#ifndef CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H
#define CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H


#include <cstdlib>


namespace cell_routines {

    namespace mpi_binary_naive {
        /**
         * Process upper triangles of squares of size sizeof(Input)*8 that lies on specific antidiagonal
         */
        template<class Input>
        inline void
        loop_upper_half_binary(int lower_bound, int upper_bound, int shift, int l_edge, int t_edge,
                               Input active_bits, Input *l_strands, Input *t_strands, Input *a_reverse, Input *b) {

            for (int j = lower_bound; j < upper_bound; ++j) {

                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input l_strand_cap = l_strand >> shift;
                Input symbol_a = a_reverse[l_edge + j];
                Input symbol_b = b[t_edge + j];
                Input cond = active_bits & (((~(symbol_a >> shift)) ^ symbol_b) | (((~(l_strand_cap)) & t_strand)));
                Input inv_cond = ~cond;

                t_strands[t_edge + j] = (inv_cond & t_strand) | (cond & l_strand_cap);
                t_strand = t_strand << shift;

                cond = cond << shift;
                inv_cond = ~cond;
                l_strands[l_edge + j] = (inv_cond & l_strand) | (cond & t_strand);
            }
        }

        /**
         * Same as loop_upper_half_binary but without offset requirement
         */
        template<class Input>
        inline void loop_center_half_binary(int lower_bound, int upper_bound, int l_edge, int t_edge,
                                            Input *l_strands, Input *t_strands, Input *a_reverse, Input *b) {
            for (int j = lower_bound; j < upper_bound; ++j) {
                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input cond = ((~(a_reverse[l_edge + j] ^ b[t_edge + j])) | ((~l_strand) & t_strand));
                Input rev_combing_cond = ~cond;

                l_strands[l_edge + j] = (rev_combing_cond & l_strand) | (cond & t_strand);
                t_strands[t_edge + j] = (rev_combing_cond & t_strand) | (cond & l_strand);
            }
        }

        template<class Input>
        inline void loop_lower_half_binary(int lower_bound, int upper_bound, int shift, int l_edge, int t_edge,
                                           Input active_bits, Input *l_strands, Input *t_strands, Input *a_reverse,
                                           Input *b) {

            for (int j = lower_bound; j < upper_bound; ++j) {

                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input l_strand_cap = l_strand << (shift + 1);
                Input symbol_a = a_reverse[l_edge + j];
                Input symbol_b = b[t_edge + j];
                Input cond =
                        active_bits & (((~(symbol_a << (shift + 1))) ^ symbol_b) | (((~(l_strand_cap)) & t_strand)));
                Input inv_cond = ~cond;

                t_strands[t_edge + j] = (inv_cond & t_strand) | (cond & l_strand_cap);

                t_strand = t_strand >> (shift + 1);

                cond = cond >> (shift + 1);
                inv_cond = ~cond;
                l_strands[l_edge + j] = (inv_cond & l_strand) | (cond & t_strand);
            }
        }

        /**
 * Process upper triangles of squares of size sizeof(Input)*8 that lies on specific antidiagonal
 */
        template<class Input>
        inline void
        loop_upper_half_binary_mpi(int lower_bound, int upper_bound, int shift, int l_edge, int t_edge,
                               Input active_bits, Input *l_strands, Input *t_strands, Input *a_reverse, Input *b) {

#pragma omp  for simd schedule(static) aligned(l_strands, t_strands, a_reverse, b:sizeof(Input)*8)
            for (int j = lower_bound; j < upper_bound; ++j) {

                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input l_strand_cap = l_strand >> shift;
                Input symbol_a = a_reverse[l_edge + j];
                Input symbol_b = b[t_edge + j];
                Input cond = active_bits & (((~(symbol_a >> shift)) ^ symbol_b) | (((~(l_strand_cap)) & t_strand)));
                Input inv_cond = ~cond;

                t_strands[t_edge + j] = (inv_cond & t_strand) | (cond & l_strand_cap);
                t_strand = t_strand << shift;

                cond = cond << shift;
                inv_cond = ~cond;
                l_strands[l_edge + j] = (inv_cond & l_strand) | (cond & t_strand);
            }
        }

        /**
         * Same as loop_upper_half_binary but without offset requirement
         */
        template<class Input>
        inline void loop_center_half_binary_mpi(int lower_bound, int upper_bound, int l_edge, int t_edge,
                                            Input *l_strands, Input *t_strands, Input *a_reverse, Input *b) {
#pragma omp  for simd schedule(static) aligned(l_strands, t_strands, a_reverse, b:sizeof(Input)*8)
            for (int j = lower_bound; j < upper_bound; ++j) {
                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input cond = ((~(a_reverse[l_edge + j] ^ b[t_edge + j])) | ((~l_strand) & t_strand));
                Input rev_combing_cond = ~cond;

                l_strands[l_edge + j] = (rev_combing_cond & l_strand) | (cond & t_strand);
                t_strands[t_edge + j] = (rev_combing_cond & t_strand) | (cond & l_strand);
            }
        }

        template<class Input>
        inline void loop_lower_half_binary_mpi(int lower_bound, int upper_bound, int shift, int l_edge, int t_edge,
                                           Input active_bits, Input *l_strands, Input *t_strands, Input *a_reverse,
                                           Input *b) {

#pragma omp  for simd schedule(static) aligned(l_strands, t_strands, a_reverse, b:sizeof(Input)*8)
            for (int j = lower_bound; j < upper_bound; ++j) {

                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input l_strand_cap = l_strand << (shift + 1);
                Input symbol_a = a_reverse[l_edge + j];
                Input symbol_b = b[t_edge + j];
                Input cond = active_bits & (((~(symbol_a << (shift + 1))) ^ symbol_b) | (((~(l_strand_cap)) & t_strand)));
                Input inv_cond = ~cond;

                t_strands[t_edge + j] = (inv_cond & t_strand) | (cond & l_strand_cap);

                t_strand = t_strand >> (shift + 1);

                cond = cond >> (shift + 1);
                inv_cond = ~cond;
                l_strands[l_edge + j] = (inv_cond & l_strand) | (cond & t_strand);
            }
        }
    }

    namespace mpi_binary_smart {}


    namespace mpi_4_size {

    }

    // generic
}

namespace prefix_lcs_via_semi_local {

    namespace  binary {
        using namespace cell_routines::mpi_binary_naive;

        template<class Input>
        int prefix_lcs_via_braid_bits_2symbol(Input *a_reverse, int a_size, Input *b, int b_size, int a_total_symbols, int threads_num = 1) {


            Input *l_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
            Input *t_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

            auto m = a_size, n = b_size;

            int dis_braid = 0;
            auto num_diag = m + n - 1;
            auto total_same_length_diag = num_diag - (m) - (m - 1);

            #pragma omp parallel num_threads(threads_num)  default(none) shared(l_strands, t_strands, a_reverse, b, m, n, dis_braid, total_same_length_diag)
            {

                Input mask;


                #pragma omp  for simd schedule(static) aligned(l_strands:sizeof(Input)*8)
                for (int k = 0; k < m; ++k) l_strands[k] = ~Input(0);

                #pragma omp  for simd schedule(static) aligned(t_strands:sizeof(Input)*8)
                for (int k = 0; k < n; ++k) t_strands[k] = Input(0);

                auto upper_bound = (sizeof(Input) * 8) - 1;

                //process first triangle in 0,0 cube
                mask = Input(1);
                Input mask_r = Input(1) << (sizeof(Input) * 8 - 1);
                int rev_counter = (sizeof(Input) * 8 - 1);

                //PHASE 0:Process first triangle
                #pragma omp single
                {
                    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter--) {
                        loop_upper_half_binary(0, 1, rev_counter, m - 1,0,mask,l_strands,t_strands,a_reverse,b);
                        mask = (mask << 1) | Input(1);
                    }
                }

                //PHASE 1: Process diagonal till fill big left triangle,
                for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {
                    //to process current
                    rev_counter = (sizeof(Input) * 8 - 1);
                    mask = Input(1);
                    //to process previous
                    Input mask_prev = ~static_cast<Input>(0);

                    //process previous size/2 - 1 cubes and current size/2 -1  cubes
                    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
                        //update mask of prev move
                        mask_prev <<= 1;

                        loop_upper_half_binary_mpi(0,cur_diag_cube_len + 1, rev_counter, m-1-cur_diag_cube_len,0,mask,l_strands,t_strands,a_reverse,b);

                        loop_lower_half_binary_mpi(0,cur_diag_cube_len,inside_diag_num,m-1-cur_diag_cube_len+1,0,mask_prev,l_strands,t_strands,a_reverse,b);


                        //update mask of current move
                        mask = (mask << 1) | Input(1);
                    }

                    loop_center_half_binary_mpi(0,cur_diag_cube_len+1,m-1-cur_diag_cube_len,0,l_strands,t_strands,a_reverse,b);

                }

                //PHASE 2
                for (int k = 0; k < total_same_length_diag; ++k) {

                    //to process current
                    rev_counter = (sizeof(Input) * 8 - 1);
                    mask = Input(1);

                    //to process previous
                    Input mask_prev = ~Input(0);

                    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
                        //update mask of prev move
                        mask_prev <<= 1;

                        loop_upper_half_binary_mpi(0, m, rev_counter, 0, k + 1, mask, l_strands,t_strands, a_reverse, b);

                        loop_lower_half_binary_mpi(0, m, inside_diag_num, 0, k, mask_prev, l_strands, t_strands, a_reverse, b);


                        //update mask of current move
                        mask = (mask << 1) | Input(1);
                    }

                    loop_center_half_binary_mpi(0, m, 0, k + 1, l_strands, t_strands, a_reverse, b);
                }


                auto start_j = total_same_length_diag + 1;
                for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

                    //to process current
                    rev_counter = (sizeof(Input) * 8 - 1);
                    mask = Input(1);

                    //to process previous
                    Input mask_prev = ~Input(0);

                    //process previous size/2 - 1 cubes and current size/2 -1  cubes
                    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
                        //update mask of prev move
                        mask_prev <<= 1;

                        loop_upper_half_binary_mpi(0, cur_diag_cube_len, rev_counter, 0, start_j, mask, l_strands, t_strands, a_reverse, b);

                        loop_lower_half_binary_mpi(0, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, mask_prev, l_strands, t_strands, a_reverse, b);
                        //update mask of current move
                        mask = (mask << 1) | Input(1);

                    }

                    loop_center_half_binary_mpi(0, cur_diag_cube_len, 0, start_j, l_strands, t_strands, a_reverse, b);
                }


                //process last triangle in position  m-1, n-1  cube
                mask = ~Input(0);
                mask_r = mask;
                #pragma omp single
                {

                    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
                        mask = mask << 1;
                        loop_lower_half_binary(0, 1, inside_diag_num, 0, n - 1, mask, l_strands, t_strands, a_reverse, b);
                    }

                }

                #pragma omp for  simd schedule(static) reduction(+:dis_braid)  aligned(l_strands:sizeof(Input)*8)
                for (int i1 = 0; i1 < m; ++i1) {
                    //  Brian Kernighan’s Algorithm
                    int counter = 0;
                    Input number = l_strands[i1];
                    //  LogNumber
                    while (number) {
                        number &= (number - 1);
                        counter++;
                    }
                    dis_braid += counter;
                }
            }

            return a_total_symbols - dis_braid;
        }
    }
}

#endif //CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H
