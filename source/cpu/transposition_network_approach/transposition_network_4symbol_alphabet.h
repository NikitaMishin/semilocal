//
// Created by nikita on 30.07.2020.
//

#ifndef CPU_TRANSPOSITION_NETWORK_4SYMBOL_ALPHABET_H
#define CPU_TRANSPOSITION_NETWORK_4SYMBOL_ALPHABET_H

#include <vector>
#include <cmath>




//2 bit per symbol
// a<=b !
//assume both multiple by sizeof(Input)*8
template<class Input>
int prefix_lcs_via_braid_bits_4symbol_mpi(Input *a_reverse, int a_size, int a_total_symbols,
                                          Input *b, int b_size, int b_total_symbols,
                                          int threads_num = 1) {


    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);

    Input braid_ones = Input(1);
    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << shift * 2);
    }


#pragma omp parallel num_threads(threads_num)  default(none) shared(bitset_left_strand_map, bitset_top_strand_map, a_reverse, b, m, n, dis_braid, total_same_length_diag, braid_ones)
    {

        int left_edge, top_edge;
        Input mask;


        #pragma omp  for simd schedule(static) nowait aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < m; ++k) {
            bitset_left_strand_map[k] = braid_ones;
        }

        #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map:sizeof(Input)*8)
        for (int k = 0; k < n; ++k) {
            bitset_top_strand_map[k] = Input(0);
        }

        auto upper_bound = (sizeof(Input) * 8 / 2) - 1;
        //process first triangle in 0,0 cube
        mask = Input(1);
        Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
        int rev_counter = (sizeof(Input) * 8 - 2);

        //PHASE 0:Process first triangle
        #pragma omp single
        {
            for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {

                Input left_strand = bitset_left_strand_map[m - 1];
                Input top_strand = bitset_top_strand_map[0];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[m - 1];
                Input symbol_b = b[0];
                Input symbols = ~((symbol_a >> rev_counter) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition =
                        mask & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[0] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);

                    top_strand = top_strand << rev_counter;
                    symbols = (~(symbol_a)) ^ (symbol_b << rev_counter);
                    symbols &= (symbols >> 1) & braid_ones;

                    combing_condition =
                            mask_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[m - 1] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);

            }
        }



        //PHASE 1: Process diagonal till fill big left triangle,
        //
        for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

            left_edge = m - 1 - cur_diag_cube_len;
            top_edge = 0;
            //to process current
            rev_counter = (sizeof(Input) * 8 - 2);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = braid_ones;
            Input mask_prev_r = braid_ones;


            //process previous size/2 - 1 cubes and current size/2 -1  cubes
            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
                //update mask of prev move
                mask_prev <<= 2;
                mask_prev_r >>= 2;

                #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand >> rev_counter;
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[j];

                    Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                    symbols &= (symbols >> 1) & braid_ones;

                    Input combing_condition =
                            mask & (symbols | (((~(left_cap)) & top_strand)));
                    Input rev_combing_cond = combing_condition ^braid_ones;


                    if (combing_condition) {
                        bitset_top_strand_map[top_edge + j] =
                                (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                        top_strand = top_strand << rev_counter;
                        symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                        symbols &= (symbols >> 1) & braid_ones;
                        combing_condition =
                                mask_r & (symbols | ((~(left_strand) & top_strand)));
                        rev_combing_cond = combing_condition ^ braid_ones;

                        bitset_left_strand_map[left_edge + j] =
                                (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    }

                }


                left_edge++;
                // parallel process previous
                #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len; ++j) {

                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[j];
                    Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                    symbols &= (symbols >> 1) & braid_ones;

                    Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                    Input rev_combing_cond = combing_condition ^braid_ones;

                    if (combing_condition) {
                        bitset_top_strand_map[top_edge + j] =
                                (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                        top_strand = top_strand >> (2 * (inside_diag_num + 1));
                        symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                        symbols &= (symbols >> 1) & braid_ones;
                        combing_condition =
                                mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                        rev_combing_cond = combing_condition ^ braid_ones;

                        bitset_left_strand_map[left_edge + j] =
                                (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    }
                }
                left_edge--;


                //update mask of current move
                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);
            }

//        process last center mask is always all ones
            // parallel process cur
            #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input symbols = (~(a_reverse[left_edge + j] ^ b[j]));
                symbols &= (symbols >> 1) & braid_ones;
                Input combing_condition = (symbols | ((~left_strand) & top_strand));
                Input rev_combing_cond = combing_condition ^braid_ones;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_strand);
                }
            }
        }



        //PHASE 2
        for (int k = 0; k < total_same_length_diag; ++k) {
            left_edge = 0;
            top_edge = k + 1;


            //to process current
            rev_counter = (sizeof(Input) * 8 - 2);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = braid_ones;
            Input mask_prev_r = braid_ones;

            //process previous size/2 - 1 cubes and current size/2 -1  cubes

            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
                //update mask of prev move
                mask_prev <<= 2;
                mask_prev_r >>= 2;

                //parallel  process current
                #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < m; ++j) {

                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand >> rev_counter;
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[top_edge + j];
                    Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                    symbols &= (symbols >> 1) & braid_ones;
                    Input combing_condition =
                            mask & (symbols | (((~(left_cap)) & top_strand)));
                    Input rev_combing_cond = combing_condition ^braid_ones;

                    if (combing_condition) {
                        bitset_top_strand_map[top_edge + j] =
                                (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                        top_strand = top_strand << rev_counter;
                        symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                        symbols &= (symbols >> 1) & braid_ones;
                        combing_condition =
                                mask_r & (symbols | ((~(left_strand) & top_strand)));
                        rev_combing_cond = combing_condition ^ braid_ones;

                        bitset_left_strand_map[left_edge + j] =
                                (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    }
                }


                top_edge--;
                // parallel process previous
                #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < m; ++j) {

                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[top_edge + j];
                    Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                    symbols &= (symbols >> 1) & braid_ones;

                    Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                    Input rev_combing_cond = combing_condition ^braid_ones;

                    if (combing_condition) {
                        bitset_top_strand_map[top_edge + j] =
                                (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                        top_strand = top_strand >> (2 * (inside_diag_num + 1));
                        symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                        symbols &= (symbols >> 1) & braid_ones;
                        combing_condition =
                                mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                        rev_combing_cond = combing_condition ^ braid_ones;

                        bitset_left_strand_map[left_edge + j] =
                                (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    }

                }
                top_edge++;


                //update mask of current move
                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);

            }




            //proccess center
            #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
            for (int j = 0; j < m; ++j) {
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input symbols = (~(a_reverse[left_edge + j] ^ b[top_edge + j]));
                symbols &= (symbols >> 1) & braid_ones;
                Input combing_condition = (symbols | ((~left_strand) & top_strand));
                Input rev_combing_cond = combing_condition ^braid_ones;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_strand);
                }


            }

        }



        //PHASE 3
        auto start_j = total_same_length_diag + 1;
        for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {
            left_edge = 0;
            top_edge = start_j;

            //to process current
            rev_counter = (sizeof(Input) * 8 - 2);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = braid_ones;
            Input mask_prev_r = braid_ones;

            //process previous size/2 - 1 cubes and current size/2 -1  cubes
            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
                //update mask of prev move
                mask_prev <<= 2;
                mask_prev_r >>= 2;

                //parallel  process current
                #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len; ++j) {
                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand >> rev_counter;
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[top_edge + j];
                    Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                    symbols &= (symbols >> 1) & braid_ones;
                    Input combing_condition =
                            mask & (symbols | (((~(left_cap)) & top_strand)));
                    Input rev_combing_cond = combing_condition ^braid_ones;

                    if (combing_condition) {
                        bitset_top_strand_map[top_edge + j] =
                                (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                        top_strand = top_strand << rev_counter;
                        symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                        symbols &= (symbols >> 1) & braid_ones;
                        combing_condition =
                                mask_r & (symbols | ((~(left_strand) & top_strand)));
                        rev_combing_cond = combing_condition ^ braid_ones;

                        bitset_left_strand_map[left_edge + j] =
                                (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    }

                }


                top_edge--;
                // parallel process previous
                #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len + 1; ++j) {


                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[top_edge + j];
                    Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                    symbols &= (symbols >> 1) & braid_ones;

                    Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                    Input rev_combing_cond = combing_condition ^braid_ones;

                    if (combing_condition) {
                        bitset_top_strand_map[top_edge + j] =
                                (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                        top_strand = top_strand >> (2 * (inside_diag_num + 1));
                        symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                        symbols &= (symbols >> 1) & braid_ones;
                        combing_condition =
                                mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                        rev_combing_cond = combing_condition ^ braid_ones;

                        bitset_left_strand_map[left_edge + j] =
                                (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    }


                }
                top_edge++;


                //update mask of current move
                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);

            }




            //proccess center
            #pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
            for (int j = 0; j < cur_diag_cube_len; ++j) {
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input symbols = (~(a_reverse[left_edge + j] ^ b[top_edge + j]));
                symbols &= (symbols >> 1) & braid_ones;
                Input combing_condition = (symbols | ((~left_strand) & top_strand));
                Input rev_combing_cond = combing_condition ^braid_ones;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_strand);
                }
            }

        }


        //process last triangle in position  m-1, n-1  cube
        mask = braid_ones;
        mask_r = mask;
        #pragma omp single
        {

            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
                mask = mask << 2;
                mask_r = mask_r >> 2;

                Input left_strand = bitset_left_strand_map[0];
                Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                Input symbol_a = a_reverse[0];
                Input symbol_b = b[n - 1];
                Input symbols = ~((symbol_a << (inside_diag_num * 2 + 2)) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input top_strand = bitset_top_strand_map[n - 1];
                Input combing_condition =
                        mask & (symbols | ((~left_cap) & top_strand));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[n - 1] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);

                    top_strand = top_strand >> (2 * (inside_diag_num + 1));
                    symbols = ~(symbol_a ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_r &
                            (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;
                    bitset_left_strand_map[0] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

            }
        }

        #pragma omp for  simd schedule(static) reduction(+:dis_braid)  aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
        for (int i1 = 0; i1 < m; ++i1) {
            //  Brian Kernighan’s Algorithm
            int counter = 0;
            Input number = bitset_left_strand_map[i1];
            //  LogNumber
            while (number) {
                number &= (number - 1);
                counter++;
            }
            dis_braid += counter;
        }
    }

    free(bitset_left_strand_map);
    free(bitset_top_strand_map);

    return a_total_symbols - dis_braid;

}


//2 bit per symbol
// a<=b !
//assume both multiple by sizeof(Input)*8
template<class Input>
int prefix_lcs_via_braid_bits_4symbol(Input *a_reverse, int a_size, int a_total_symbols,
                                          Input *b, int b_size, int b_total_symbols) {

    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);

    Input braid_ones = Input(1);
    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << shift * 2);
    }


    int left_edge, top_edge;
    Input mask;

    for (int k = 0; k < m; ++k) {
        bitset_left_strand_map[k] = braid_ones;
    }

    for (int k = 0; k < n; ++k) {
        bitset_top_strand_map[k] = Input(0);
    }

    auto upper_bound = (sizeof(Input) * 8 / 2) - 1;
    //process first triangle in 0,0 cube
    mask = Input(1);
    Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
    int rev_counter = (sizeof(Input) * 8 - 2);

    //PHASE 0:Process first triangle
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {

        Input left_strand = bitset_left_strand_map[m - 1];
        Input top_strand = bitset_top_strand_map[0];
        Input left_cap = left_strand >> rev_counter;
        Input symbol_a = a_reverse[m - 1];
        Input symbol_b = b[0];
        Input symbols = ~((symbol_a >> rev_counter) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        Input combing_condition =
                mask & (symbols | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = combing_condition ^braid_ones;

        if (combing_condition) {
            bitset_top_strand_map[0] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            top_strand = top_strand << rev_counter;
            symbols = (~(symbol_a)) ^ (symbol_b << rev_counter);
            symbols &= (symbols >> 1) & braid_ones;

            combing_condition =
                    mask_r & (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            bitset_left_strand_map[m - 1] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);

    }



    //PHASE 1: Process diagonal till fill big left triangle,
    //
    for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

        left_edge = m - 1 - cur_diag_cube_len;
        top_edge = 0;
        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        //process previous size/2 - 1 cubes and current size/2 -1  cubes
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[j];

                Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition =
                        mask & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;


                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand << rev_counter;
                    symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

            }


            left_edge++;
            // parallel process previous
            for (int j = 0; j < cur_diag_cube_len; ++j) {

                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[j];
                Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand >> (2 * (inside_diag_num + 1));
                    symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }
            }
            left_edge--;


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

//        process last center mask is always all ones
        // parallel process cur
        for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input symbols = (~(a_reverse[left_edge + j] ^ b[j]));
            symbols &= (symbols >> 1) & braid_ones;
            Input combing_condition = (symbols | ((~left_strand) & top_strand));
            Input rev_combing_cond = combing_condition ^braid_ones;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            }
        }
    }



    //PHASE 2
    for (int k = 0; k < total_same_length_diag; ++k) {
        left_edge = 0;
        top_edge = k + 1;


        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;

        //process previous size/2 - 1 cubes and current size/2 -1  cubes

        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            //parallel  process current
            for (int j = 0; j < m; ++j) {

                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;
                Input combing_condition =
                        mask & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand << rev_counter;
                    symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }
            }


            top_edge--;
            // parallel process previous
            for (int j = 0; j < m; ++j) {

                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand >> (2 * (inside_diag_num + 1));
                    symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

            }
            top_edge++;


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);

        }




        //proccess center
        for (int j = 0; j < m; ++j) {
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input symbols = (~(a_reverse[left_edge + j] ^ b[top_edge + j]));
            symbols &= (symbols >> 1) & braid_ones;
            Input combing_condition = (symbols | ((~left_strand) & top_strand));
            Input rev_combing_cond = combing_condition ^braid_ones;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            }


        }

    }



    //PHASE 3
    auto start_j = total_same_length_diag + 1;
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {
        left_edge = 0;
        top_edge = start_j;

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;

        //process previous size/2 - 1 cubes and current size/2 -1  cubes
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            //parallel  process current
            for (int j = 0; j < cur_diag_cube_len; ++j) {
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;
                Input combing_condition =
                        mask & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand << rev_counter;
                    symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

            }


            top_edge--;
            // parallel process previous
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {


                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand >> (2 * (inside_diag_num + 1));
                    symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }


            }
            top_edge++;


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);

        }




        //proccess center
        for (int j = 0; j < cur_diag_cube_len; ++j) {
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input symbols = (~(a_reverse[left_edge + j] ^ b[top_edge + j]));
            symbols &= (symbols >> 1) & braid_ones;
            Input combing_condition = (symbols | ((~left_strand) & top_strand));
            Input rev_combing_cond = combing_condition ^braid_ones;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            }
        }

    }


    //process last triangle in position  m-1, n-1  cube
    mask = braid_ones;
    mask_r = mask;

    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        Input left_strand = bitset_left_strand_map[0];
        Input left_cap = left_strand << (2 * (inside_diag_num + 1));
        Input symbol_a = a_reverse[0];
        Input symbol_b = b[n - 1];
        Input symbols = ~((symbol_a << (inside_diag_num * 2 + 2)) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        Input top_strand = bitset_top_strand_map[n - 1];
        Input combing_condition =
                mask & (symbols | ((~left_cap) & top_strand));
        Input rev_combing_cond = combing_condition ^braid_ones;

        if (combing_condition) {
            bitset_top_strand_map[n - 1] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            top_strand = top_strand >> (2 * (inside_diag_num + 1));
            symbols = ~(symbol_a ^ (symbol_b >> (2 * (inside_diag_num + 1))));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition =
                    mask_r &
                    (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;
            bitset_left_strand_map[0] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }

    }
    for (int i1 = 0; i1 < m; ++i1) {
        //  Brian Kernighan’s Algorithm
        int counter = 0;
        Input number = bitset_left_strand_map[i1];
        //  LogNumber
        while (number) {
            number &= (number - 1);
            counter++;
        }
        dis_braid += counter;
    }

    free(bitset_left_strand_map);
    free(bitset_top_strand_map);
    
    return a_total_symbols - dis_braid;

}






//2 bit per symbol
// a<=b !
//assume both multiple by sizeof(Input)*8
template<class Input>
int prefix_lcs_via_braid_bits_4symbol_finale(Input *a_reverse, int a_size, int a_total_symbols,
                                      Input *b, int b_size, int b_total_symbols) {

    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);

    Input braid_ones = Input(1);
    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << shift * 2);
    }

    int active_symbols_a_active = a_total_symbols % (sizeof(Input)*8/2);
    int active_symbols_b_active = b_total_symbols % (sizeof(Input)*8/2);
    Input l_active_mask = active_symbols_a_active ? Input(0) : braid_ones;
    Input r_active_mask = active_symbols_b_active ? Input(0) : braid_ones;




    for (int i = 0; i <active_symbols_a_active ; ++i) {
        l_active_mask|= (Input(1)<<((sizeof(Input) * 8-2)-i * 2));
    }

    for (int i = 0; i <active_symbols_b_active ; ++i) {
        r_active_mask |= (Input(1) << (2 * i));
    }

    int left_edge, top_edge;
    Input mask;

    for (int k = 0; k < m; ++k) {
        bitset_left_strand_map[k] = braid_ones;
    }


    for (int k = 0; k < n; ++k) {
        bitset_top_strand_map[k] = Input(0);
    }

    auto upper_bound = (sizeof(Input) * 8 / 2) - 1;
    //process first triangle in 0,0 cube
    mask = Input(1);
    Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
    int rev_counter = (sizeof(Input) * 8 - 2);


    //PHASE 0:Process first triangle
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {

        Input left_strand = bitset_left_strand_map[m - 1];
        Input top_strand = bitset_top_strand_map[0];
        Input left_cap = left_strand >> rev_counter;
        Input symbol_a = a_reverse[m - 1];
        Input symbol_b = b[0];
        Input symbols = ~((symbol_a >> rev_counter) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        Input combing_condition =
                mask & (symbols | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = combing_condition ^braid_ones;

        if (combing_condition) {
            bitset_top_strand_map[0] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            top_strand = top_strand << rev_counter;
            symbols = (~(symbol_a)) ^ (symbol_b << rev_counter);
            symbols &= (symbols >> 1) & braid_ones;

            combing_condition =
                    mask_r & (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            bitset_left_strand_map[m - 1] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);

    }



    //PHASE 1: Process diagonal till fill big left triangle,
    //
    for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

        left_edge = m - 1 - cur_diag_cube_len;
        top_edge = 0;
        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        //process previous size/2 - 1 cubes and current size/2 -1  cubes
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[j];

                Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition =
                        mask & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;


                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand << rev_counter;
                    symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

            }


            left_edge++;
            // parallel process previous
            for (int j = 0; j < cur_diag_cube_len; ++j) {

                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[j];
                Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand >> (2 * (inside_diag_num + 1));
                    symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }
            }
            left_edge--;


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

//        process last center mask is always all ones
        // parallel process cur
        for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input symbols = (~(a_reverse[left_edge + j] ^ b[j]));
            symbols &= (symbols >> 1) & braid_ones;
            Input combing_condition = (symbols | ((~left_strand) & top_strand));
            Input rev_combing_cond = combing_condition ^braid_ones;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            }
        }
    }



    //PHASE 2
    for (int k = 0; k < total_same_length_diag; ++k) {
        left_edge = 0;
        top_edge = k + 1;


        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;

        //process previous size/2 - 1 cubes and current size/2 -1  cubes

        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            //parallel  process current
            for (int j = 0; j < m; ++j) {

                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;
                Input combing_condition =
                        mask & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand << rev_counter;
                    symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }
            }


            top_edge--;
            // parallel process previous
            for (int j = 0; j < m; ++j) {

                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand >> (2 * (inside_diag_num + 1));
                    symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

            }
            top_edge++;


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);

        }




        //proccess center
        for (int j = 0; j < m; ++j) {
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input symbols = (~(a_reverse[left_edge + j] ^ b[top_edge + j]));
            symbols &= (symbols >> 1) & braid_ones;
            Input combing_condition = (symbols | ((~left_strand) & top_strand));
            Input rev_combing_cond = combing_condition ^braid_ones;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            }


        }

    }



    //PHASE 3
    auto start_j = total_same_length_diag + 1;
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {
        left_edge = 0;
        top_edge = start_j;

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;

        //process previous size/2 - 1 cubes and current size/2 -1  cubes
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            //parallel  process current
            for (int j = 0; j < cur_diag_cube_len; ++j) {
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;
                Input combing_condition =
                        mask & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand << rev_counter;
                    symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

            }


            top_edge--;
            // parallel process previous
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {


                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand << (2 * (inside_diag_num + 1));
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[top_edge + j];
                Input symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
                symbols &= (symbols >> 1) & braid_ones;

                Input combing_condition = mask_prev & (symbols | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = combing_condition ^braid_ones;

                if (combing_condition) {
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);
                    top_strand = top_strand >> (2 * (inside_diag_num + 1));
                    symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
                    symbols &= (symbols >> 1) & braid_ones;
                    combing_condition =
                            mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
                    rev_combing_cond = combing_condition ^ braid_ones;

                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }


            }
            top_edge++;


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);

        }




        //proccess center
        for (int j = 0; j < cur_diag_cube_len; ++j) {
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input symbols = (~(a_reverse[left_edge + j] ^ b[top_edge + j]));
            symbols &= (symbols >> 1) & braid_ones;
            Input combing_condition = (symbols | ((~left_strand) & top_strand));
            Input rev_combing_cond = combing_condition ^braid_ones;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            }
        }

    }


    //process last triangle in position  m-1, n-1  cube
    mask = braid_ones;
    mask_r = mask;

    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        Input left_strand = bitset_left_strand_map[0];
        Input left_cap = left_strand << (2 * (inside_diag_num + 1));
        Input symbol_a = a_reverse[0];
        Input symbol_b = b[n - 1];
        Input symbols = ~((symbol_a << (inside_diag_num * 2 + 2)) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        Input top_strand = bitset_top_strand_map[n - 1];
        Input combing_condition =
                mask & (symbols | ((~left_cap) & top_strand));
        Input rev_combing_cond = combing_condition ^braid_ones;

        if (combing_condition) {
            bitset_top_strand_map[n - 1] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            top_strand = top_strand >> (2 * (inside_diag_num + 1));
            symbols = ~(symbol_a ^ (symbol_b >> (2 * (inside_diag_num + 1))));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition =
                    mask_r &
                    (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;
            bitset_left_strand_map[0] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }

    }
    for (int i1 = 0; i1 < m; ++i1) {
        //  Brian Kernighan’s Algorithm
        int counter = 0;
        Input number = bitset_left_strand_map[i1];
        //  LogNumber
        while (number) {
            number &= (number - 1);
            counter++;
        }
        dis_braid += counter;
    }

    free(bitset_left_strand_map);
    free(bitset_top_strand_map);

    return a_total_symbols - dis_braid;

}


template<class Input>
int prefix_lcs_via_braid_4symbol_one_one_size(Input a_reverse, int a_total_symbols, Input b, int b_total_symbols) {

    int dis_braid = 0;
    Input braid_ones = Input(1);
    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << (shift * 2));
    }

    int active_symbols_a_active = a_total_symbols % (sizeof(Input)*8/2);
    int active_symbols_b_active = b_total_symbols % (sizeof(Input)*8/2);
    Input l_active_mask = active_symbols_a_active ? Input(0) : braid_ones;
    Input r_active_mask = active_symbols_b_active ? Input(0) : braid_ones;
    for (int i = 0; i <active_symbols_a_active ; ++i) {
        l_active_mask|= (Input(1)<<((sizeof(Input) * 8-2)-i * 2));
    }
    for (int i = 0; i <active_symbols_b_active ; ++i) {
        r_active_mask |= (Input(1) << (2 * i));
    }

    a_reverse <<= (sizeof(Input)*8 - 2*active_symbols_a_active ); //since we have a_reverse =  00000 last,.....


    int upper_bound = (sizeof(Input) * 8 / 2) - 1;
    Input left_strand = l_active_mask;
    Input top_strand = Input(0);
    Input old_top;
    Input mask = Input(1);
    Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
    int rev_counter = (sizeof(Input) * 8 - 2);

    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {
        Input left_cap = left_strand >> rev_counter;
        Input symbols = ~((a_reverse >> rev_counter) ^ b);
        symbols &= (symbols >> 1) & braid_ones;

        Input combing_condition = (r_active_mask) & (l_active_mask>>rev_counter) & mask & (symbols | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = combing_condition ^ braid_ones;

        if (combing_condition) {
            old_top = top_strand << rev_counter;
            top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            symbols = (~(a_reverse)) ^ (b << rev_counter);
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = l_active_mask & (r_active_mask<<rev_counter) & mask_r & (symbols | ((~(left_strand) & old_top)));
            rev_combing_cond = combing_condition ^ braid_ones;
            left_strand = (rev_combing_cond & left_strand) | (combing_condition & old_top);
        }

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);
    }

    mask = braid_ones;
    mask_r = mask;
    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        Input left_cap = left_strand << (2 * (inside_diag_num + 1));
        Input symbols = ~((a_reverse << (inside_diag_num * 2 + 2)) ^ b);
        symbols &= (symbols >> 1) & braid_ones;

        Input combing_condition = (l_active_mask << (2 * (inside_diag_num + 1))) & r_active_mask & mask & (symbols | ((~left_cap) & top_strand));
        Input rev_combing_cond = combing_condition ^ braid_ones;

        if (combing_condition) {
            old_top = top_strand >> (2 * (inside_diag_num + 1));
            top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            symbols = ~(a_reverse ^ (b >> (2 * (inside_diag_num + 1))));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = l_active_mask &  (r_active_mask>>(2 * (inside_diag_num + 1))) & mask_r & (symbols | ((~(left_strand) & old_top)));
            rev_combing_cond = combing_condition ^ braid_ones;
            left_strand = (rev_combing_cond & left_strand) | (combing_condition & old_top);
        }
    }




    int counter = 0;
    while (left_strand) {
        left_strand &= (left_strand - 1);
        counter++;
    }
    dis_braid += counter;


    return a_total_symbols - dis_braid;
}



#endif //CPU_TRANSPOSITION_NETWORK_4SYMBOL_ALPHABET_H
