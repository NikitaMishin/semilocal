//
// Created by nikita on 19.08.2020.
//

#ifndef CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H
#define CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H


#include <cstdlib>

/**
 * Contains  implementations of cell processing logic for bitwise algorithm based on semi-local iterative combing
 */
namespace cell_routines {

    namespace mpi_binary_naive {
        /**
         * Process upper triangles of squares of size sizeof(Input)*8 that lies on specific antidiagonal
         *  First version
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

    namespace mpi_binary_smart {

        template<class Input>
        inline void process_antidiag_formula1(int lower_bound, int upper_bound, int l_edge, int t_edge,
                                         Input *l_strands, Input *t_strands, Input const *a_reverse, Input const *b) {

            const int strands_per_word = sizeof(Input) * 8 - 1;

            /**
             * The idea is as follows.
             * to process some antidiagonal that have been built upon Input cells (that contains a batch of strands) we have to
             * process each such square (Input \times Input) in the antidiagonal fashion.
             * While processing each antidiagonal of such square, we need to implement the following logic:
             * in each iteration swap only those bit-strands that  active in iteration &  satisfy the combing condition.
             * This logic can be implemented by several formulas, here is presented one of it.
             * There is 20 operations inside cycle
             */
            #pragma omp   for  simd schedule(static)  aligned(l_strands, t_strands:sizeof(Input)*8) aligned(a_reverse, b:sizeof(Input)*8)
            for (int j = lower_bound; j < upper_bound; ++j) {

                Input l_strand_cap, cond, inv_cond, t_strand_shifted;
                //load phase
                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input symbol_a = a_reverse[l_edge + j];
                Input symbol_b = b[t_edge + j];

                Input mask = Input(1);

                // manual say 256 just for complete
                #pragma GCC unroll  256
                for (int shift = strands_per_word; shift > 0; shift--) {


                    l_strand_cap = l_strand >> shift;

                    cond = mask & ((~(((symbol_a >> shift)) ^ symbol_b)) | (((~(l_strand_cap)) & t_strand)));
                    inv_cond = ~cond;

                    t_strand_shifted = t_strand << shift;
                    t_strand = (inv_cond & t_strand) | (cond & l_strand_cap);

                    cond <<= shift;
                    inv_cond = ~cond;

                    l_strand = (inv_cond & l_strand) | (cond & t_strand_shifted);

                    mask = (mask << 1) | Input(1);
                }

                // center
                cond = (~(symbol_a ^ symbol_b));
                cond = (cond | ((~l_strand) & t_strand));
                inv_cond = ~cond;
                t_strand_shifted = t_strand;
                t_strand = (inv_cond & t_strand) | (cond & l_strand);
                l_strand = (inv_cond & l_strand) | (cond & t_strand_shifted);

                mask = ~Input(0);

                //lower half
                #pragma GCC unroll 256
                for (int shift = 1; shift < strands_per_word + 1; shift++) {
                    mask <<= 1;

                    l_strand_cap = l_strand << (shift);
                    cond = ~(((symbol_a << ((shift))) ^ symbol_b)); // NOT A XOR B = NOT (A XOR B)// ECONOMY

                    cond = mask & (cond | (((~(l_strand_cap)) & t_strand)));
                    inv_cond = ~cond;

                    t_strand_shifted = t_strand >> shift;
                    t_strand = (inv_cond & t_strand) | (cond & l_strand_cap);
                    cond >>= shift;
                    inv_cond = ~cond;

                    l_strand = (inv_cond & l_strand) | (cond & t_strand_shifted);
                }

                l_strands[l_edge + j] = l_strand;
                t_strands[t_edge + j] = t_strand;
            }
        }


        template<class Input>
        inline void process_antidiag_formula2(int lower_bound, int upper_bound, int l_edge, int t_edge,
                                         Input *l_strands, Input *t_strands, Input const *a_reverse, Input const *b) {

            const int strands_per_word = sizeof(Input) * 8 - 1;

            /**
             * The idea is as follows.
             * to process some antidiagonal that have been built upon Input cells (that contains a batch of strands) we have to
             * process each such square (Input \times Input) in the antidiagonal fashion.
             * While processing each antidiagonal of such square, we need to implement the following logic:
             * in each iteration swap only those bit-strands that  active in iteration &  satisfy the combing condition.
             * This logic can be implemented by several formulas, here is presented one of it.
             * There is 15 operations inside cycle
             */
            #pragma omp   for  simd schedule(static)  aligned(l_strands, t_strands:sizeof(Input)*8) aligned(a_reverse, b:sizeof(Input)*8)
            for (int j = lower_bound; j < upper_bound; ++j) {

                Input t_strand_cap;
                Input l_strand_cap, cond;

                //load
                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input symbol_a = a_reverse[l_edge + j];
                Input symbol_b = b[t_edge + j];

                Input mask = Input(1);

                // manual say 256 just for complete
                #pragma GCC unroll  256
                for (int shift = strands_per_word; shift > 0; shift--) {
                    // 15 operations inside cycle
                    // could be reduced to 14 if we store not a but ~a in memory

                    l_strand_cap = l_strand >> shift;
                    t_strand_cap = t_strand << shift;

                    cond = ~((symbol_a >> shift) ^ symbol_b);

                    t_strand = (l_strand_cap | (~mask)) & (t_strand | ( cond & mask));
                    l_strand = t_strand_cap ^ (t_strand << shift) ^ l_strand;

                    mask = (mask << 1) | Input(1);
                }

                // center, no shifts
                cond = ~((symbol_a ^ symbol_b));
                l_strand_cap = l_strand;
                t_strand_cap = t_strand;

                t_strand = (l_strand_cap | (~mask)) & (t_strand | ( cond & mask));
                l_strand = t_strand_cap ^ (t_strand) ^ l_strand;

                mask = ~Input(0);

                //lower half
                #pragma GCC unroll 256
                for (int shift = 1; shift < strands_per_word + 1; shift++) {
                    mask <<= 1;

                    l_strand_cap = l_strand << shift;
                    t_strand_cap = t_strand >> shift;
                    cond = ~(((symbol_a << (shift)) ^ symbol_b));
                    t_strand = (l_strand_cap | (~mask)) & (t_strand | ( cond & mask));
                    l_strand = t_strand_cap ^ (t_strand >> shift) ^ l_strand;
                }

                // store
                l_strands[l_edge + j] = l_strand;
                t_strands[t_edge + j] = t_strand;
            }
        }


    }


    namespace mpi_nary_size {
        template<class Input>
        inline void process_antidiagonal(int lower_bound, int upper_bound, int l_edge, int t_edge,
                                         Input *l_strands, Input *t_strands, Input const *a_reverse, Input const *b,
                                         int residue,int bits_per_strand, Input braid_ones) {
            Input single_strand = Input(1) << residue;
            int size = sizeof(Input) * 8 - residue - bits_per_strand;

            #pragma omp   for  simd schedule(static)  aligned(l_strands, t_strands:sizeof(Input)*8) aligned(a_reverse, b:sizeof(Input)*8)
            for (int j = lower_bound; j < upper_bound; ++j) {

                Input l_strand_cap, cond, t_strand_cap,eq;
                //load phase
                Input l_strand = l_strands[l_edge + j];
                Input t_strand = t_strands[t_edge + j];
                Input symbol_a = a_reverse[l_edge + j];
                Input symbol_b = b[t_edge + j];

                Input mask = single_strand;

                // manual say 256 just for complete
                #pragma GCC unroll  256
                for (int shift = size; shift > 0; shift -= bits_per_strand) {


                    l_strand_cap = l_strand >> shift;
                    t_strand_cap = t_strand << shift;

                    //reduction block
                    cond =  ~( ((symbol_a >> shift)  ) ^ symbol_b);
                    eq = cond;
                    #pragma GCC unroll  10
                    for(int i = 1; i < bits_per_strand; i++) {
                        cond &= (eq >> i);
                    }

                    t_strand = (l_strand_cap | (braid_ones ^ mask)) & (t_strand | ( cond & mask));
                    l_strand = t_strand_cap ^ (t_strand << shift) ^ l_strand;

                    mask = (mask << bits_per_strand) | single_strand;

                 }

                cond =  ~( (symbol_a) ^ symbol_b);
                eq = cond;

                #pragma GCC unroll  10
                for(int i = 1; i < bits_per_strand; i++) cond &= (eq >> i);

                l_strand_cap = l_strand;
                t_strand_cap = t_strand;


                t_strand = (l_strand_cap | (braid_ones ^ braid_ones)) & (t_strand | ( cond & braid_ones));


                l_strand = t_strand_cap ^ t_strand ^ l_strand;

                mask = braid_ones << bits_per_strand;

                //lower half
                #pragma GCC unroll 256
                for (int shift = bits_per_strand; shift <= size ; shift += bits_per_strand) {

                    //reduction block
                    cond =  ~((symbol_a << shift) ^ symbol_b);
                    eq = cond;
                    #pragma GCC unroll  10
                    for(int i = 1; i < bits_per_strand; i++) cond &= (eq >> i);

                    l_strand_cap = l_strand << shift;
                    t_strand_cap = t_strand >> shift;

                    t_strand = (l_strand_cap | (braid_ones ^ mask)) & (t_strand | ( cond & mask));
                    l_strand = (t_strand_cap ^ (t_strand >> shift) ^ l_strand);

                    mask <<= bits_per_strand;
                }

                l_strands[l_edge + j] = l_strand;
                t_strands[t_edge + j] = t_strand;
            }
        }

    }

    // generic
}

namespace prefix_lcs_via_semi_local {

    /**
     * Contains an bit-wise implementations of semi-local lcs that can compute length of LCS (llcs) of two binary string
     */
    namespace  binary {
        /**
         * This is the first non-optimized version of algorithm
         * Several condition should be satisfied:
         * 1) Since both strings are binary the elements should be compressed and stored in bits
         * 2) The first string assumed be less then the second one
         * 3) First string stored in reverse order to allow us to access elements in cell processing routine in sequential manner
         * 4) Input, that store pack of bits, should be unsigned to eliminate problem with signed shift
         * 5) Size of strings should be divisible by sizeof(Input)*8, if no, see details in paper or implementation
         * Algorithm follows the idea of iterative combing  algorithm but strands have only 0 and 1 number.
         * To see  cell processing routine see documentation of methods that used within the algorithm.
         * @tparam Input
         * @param a_reverse
         * @param a_size
         * @param b
         * @param b_size
         * @param a_total_symbols
         * @param threads_num
         * @return
         */
        template<class Input>
        int llcs_2symbol_naive_combing(Input *a_reverse, int a_size, Input *b, int b_size, int a_total_symbols, int threads_num = 1) {
            using namespace cell_routines::mpi_binary_naive;

            // also stored in the reverse order
            Input *l_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
            Input *t_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

            auto m = a_size, n = b_size;

            // total amount of strands that at the end hit right border of grid
            int dis_braid = 0;

            auto num_diag = m + n - 1;

            auto total_same_length_diag = num_diag - (m) - (m - 1);

            #pragma omp parallel num_threads(threads_num)  default(none) shared(l_strands, t_strands, a_reverse, b, m, n, dis_braid, total_same_length_diag)
            {

                Input mask;
                // Initialization step, strands that hit the left grid border all have number 1; strands that hit top grid border are 0.
                #pragma omp  for simd schedule(static) aligned(l_strands:sizeof(Input)*8)
                for (int k = 0; k < m; ++k) l_strands[k] = ~Input(0);
                #pragma omp  for simd schedule(static) aligned(t_strands:sizeof(Input)*8)
                for (int k = 0; k < n; ++k) t_strands[k] = Input(0);

                auto upper_bound = (sizeof(Input) * 8) - 1;

                //process first triangle in 0,0 cube
                mask = Input(1);
                Input mask_r = Input(1) << (sizeof(Input) * 8 - 1);
                int bits_shift = (sizeof(Input) * 8 - 1);

                //PHASE 0:Process first triangle
                #pragma omp single
                {
                    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, bits_shift--) {
                        loop_upper_half_binary(0, 1, bits_shift, m - 1, 0, mask, l_strands, t_strands, a_reverse, b);
                        mask = (mask << 1) | Input(1);
                    }
                }

                //PHASE 1: Process diagonal till fill big left triangle,
                for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {
                    //to process current
                    bits_shift = (sizeof(Input) * 8 - 1);
                    mask = Input(1);
                    //to process previous
                    Input mask_prev = ~static_cast<Input>(0);

                    //process previous size/2 - 1 cubes and current size/2 -1  cubes
                    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, bits_shift--) {
                        //update mask of prev move
                        mask_prev <<= 1;

                        loop_upper_half_binary_mpi(0,cur_diag_cube_len + 1, bits_shift, m - 1 - cur_diag_cube_len, 0, mask, l_strands, t_strands, a_reverse, b);

                        loop_lower_half_binary_mpi(0,cur_diag_cube_len,inside_diag_num,m-1-cur_diag_cube_len+1,0,mask_prev,l_strands,t_strands,a_reverse,b);


                        //update mask of current move
                        mask = (mask << 1) | Input(1);
                    }

                    loop_center_half_binary_mpi(0,cur_diag_cube_len+1,m-1-cur_diag_cube_len,0,l_strands,t_strands,a_reverse,b);

                }

                //PHASE 2
                for (int k = 0; k < total_same_length_diag; ++k) {

                    //to process current
                    bits_shift = (sizeof(Input) * 8 - 1);
                    mask = Input(1);

                    //to process previous
                    Input mask_prev = ~Input(0);

                    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, bits_shift--) {
                        //update mask of prev move
                        mask_prev <<= 1;

                        loop_upper_half_binary_mpi(0, m, bits_shift, 0, k + 1, mask, l_strands, t_strands, a_reverse, b);

                        loop_lower_half_binary_mpi(0, m, inside_diag_num, 0, k, mask_prev, l_strands, t_strands, a_reverse, b);


                        //update mask of current move
                        mask = (mask << 1) | Input(1);
                    }

                    loop_center_half_binary_mpi(0, m, 0, k + 1, l_strands, t_strands, a_reverse, b);
                }


                auto start_j = total_same_length_diag + 1;
                for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

                    //to process current
                    bits_shift = (sizeof(Input) * 8 - 1);
                    mask = Input(1);

                    //to process previous
                    Input mask_prev = ~Input(0);

                    //process previous size/2 - 1 cubes and current size/2 -1  cubes
                    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, bits_shift--) {
                        //update mask of prev move
                        mask_prev <<= 1;

                        loop_upper_half_binary_mpi(0, cur_diag_cube_len, bits_shift, 0, start_j, mask, l_strands, t_strands, a_reverse, b);

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



        template<class Input>
        int llcs_2symbol_smart_combing(Input *a_reverse, int a_size, Input *b, int b_size, int a_total_symbols, int threads_num = 1, bool formula_one = false) {
            using namespace cell_routines::mpi_binary_smart;

            // also stored in the reverse order
            Input *l_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
            Input *t_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

            auto m = a_size, n = b_size;

            // total amount of strands that at the end hit right border of grid
            int dis_braid = 0;

            auto num_diag = m + n - 1;

            auto total_same_length_diag = num_diag - (m-1) - (m - 1);

            #pragma omp parallel num_threads(threads_num)  default(none) shared(l_strands, t_strands, a_reverse, b, m, n, dis_braid, total_same_length_diag, formula_one)
            {

                // Initialization step, strands that hit the left grid border all have number 1; strands that hit top grid border are 0.
                #pragma omp  for simd schedule(static) aligned(l_strands:sizeof(Input)*8)
                for (int k = 0; k < m; ++k) l_strands[k] = ~Input(0);
                #pragma omp  for simd schedule(static) aligned(t_strands:sizeof(Input)*8)
                for (int k = 0; k < n; ++k) t_strands[k] = Input(0);

                // phase 1: process upper left triangle
                for (int diag_len = 0; diag_len < m - 1; diag_len++) {
                    formula_one?
                        process_antidiag_formula1(0,diag_len + 1,m - 1 - diag_len,0,l_strands,t_strands,a_reverse,b):
                        process_antidiag_formula2(0,diag_len + 1,m - 1 - diag_len,0,l_strands,t_strands,a_reverse,b);
                }

                // phase2: process parallelogram
                for (int k = 0; k < total_same_length_diag; k++) {
                    formula_one ?
                        process_antidiag_formula1(0, m, 0, k, l_strands, t_strands, a_reverse, b) :
                        process_antidiag_formula2(0, m, 0, k, l_strands, t_strands, a_reverse, b);
                }

                auto start_j = total_same_length_diag;

                // phase:3: lower-right triangle
                for (int diag_len = m - 1; diag_len >= 1; diag_len--) {
                    formula_one ?
                        process_antidiag_formula1(0, diag_len, 0, start_j, l_strands, t_strands, a_reverse, b) :
                        process_antidiag_formula2(0, diag_len, 0, start_j, l_strands, t_strands, a_reverse, b);
                    start_j++;
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

    namespace  nary {
        template<class Input>
        int llcs_nary_symbol_smart_combing(Input *a_reverse, int a_size, Input *b, int b_size,
                                           int a_total_symbols, int bits_per_symbol, int threads_num = 1) {
            using namespace cell_routines::mpi_nary_size;

            std::cout<<a_size<<','<<b_size<<std::endl;
            Input *l_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
            Input *t_strands = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

            auto m = a_size, n = b_size;

            // total amount of strands that at the end hit right border of grid
            int dis_braid = 0;

            auto num_diag = m + n - 1;

            auto total_same_length_diag = num_diag - (m - 1) - (m - 1);

            int residue = (sizeof(Input) * 8) % bits_per_symbol;
            Input braid_ones = 1 << residue;
            for (int i = 0; i < sizeof(Input) * 8; i += bits_per_symbol) braid_ones |= (braid_ones << i);

            #pragma omp parallel num_threads(threads_num)  default(none) shared(std::cout,residue,bits_per_symbol,braid_ones,l_strands, t_strands, a_reverse, b, m, n, dis_braid, total_same_length_diag)
            {

                // Initialization step, strands that hit the left grid border all have number 1; strands that hit top grid border are 0.
                #pragma omp  for simd schedule(static) aligned(l_strands:sizeof(Input)*8)
                for (int k = 0; k < m; ++k) l_strands[k] = braid_ones;
                #pragma omp  for simd schedule(static) aligned(t_strands:sizeof(Input)*8)
                for (int k = 0; k < n; ++k) t_strands[k] = Input(0);

                // phase 1: process upper left triangle
                for (int diag_len = 0; diag_len < m - 1; diag_len++) {
                    process_antidiagonal(0, diag_len + 1, m - 1 - diag_len, 0, l_strands, t_strands, a_reverse, b,
                                         residue, bits_per_symbol,braid_ones);
                }

                // phase2: process parallelogram
                for (int k = 0; k < total_same_length_diag; k++) {
                    process_antidiagonal(0, m, 0, k, l_strands, t_strands, a_reverse, b,residue,bits_per_symbol,braid_ones);
                }

                auto start_j = total_same_length_diag;

                // phase:3: lower-right triangle
                for (int diag_len = m - 1; diag_len >= 1; diag_len--) {
                    process_antidiagonal(0, diag_len, 0, start_j, l_strands, t_strands, a_reverse, b,
                                     residue, bits_per_symbol, braid_ones);
                    start_j++;
                }

                #pragma omp for  simd schedule(static) reduction(+:dis_braid)  aligned(l_strands:sizeof(Input)*8)
                for (int i1 = 0; i1 < m; ++i1) {
//
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
