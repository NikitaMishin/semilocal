//
// Created by nikita on 21.08.2020.
//

#ifndef CPU_TMP_H
#define CPU_TMP_H


#include <vector>
#include <cmath>

template<class Input>
inline void loop_upper_half(int lower_bound, int upper_bound, int rev_counter, int left_edge, int top_edge,
                            Input braid_ones,
                            Input mask,
                            Input mask_r,
                            Input *bitset_left_strand_map,
                            Input *bitset_top_strand_map,
                            Input *a_reverse, Input *b) {

    for (int j = lower_bound; j < upper_bound; ++j) {
        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input left_cap = left_strand >> rev_counter;
        Input symbol_a = a_reverse[left_edge + j];
        Input symbol_b = b[top_edge + j];

        Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        Input combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = combing_condition ^braid_ones;

        if (combing_condition) {
            bitset_top_strand_map[top_edge + j] = (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            top_strand = top_strand << rev_counter;
            symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = mask_r & (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }
    }
}

template<class Input>
inline void loop_center_half(int lower_bound, int upper_bound,
                             int left_edge, int top_edge,
                             Input braid_ones,
                             Input *bitset_left_strand_map,
                             Input *bitset_top_strand_map,
                             Input *a_reverse, Input *b) {

    for (int j = lower_bound; j < upper_bound; ++j) {
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

template<class Input>
inline void loop_center_half_with_mask(int lower_bound, int upper_bound,
                                       int left_edge, int top_edge,
                                       Input braid_ones, Input mask,
                                       Input *bitset_left_strand_map,
                                       Input *bitset_top_strand_map,
                                       Input *a_reverse, Input *b) {

    for (int j = lower_bound; j < upper_bound; ++j) {
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input symbols = (~(a_reverse[left_edge + j] ^ b[top_edge + j]));
        symbols &= (symbols >> 1) & braid_ones;
        Input combing_condition = mask & (symbols | ((~left_strand) & top_strand));
        Input rev_combing_cond = combing_condition ^braid_ones;
        if (combing_condition) {
            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
            bitset_top_strand_map[top_edge + j] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_strand);
        }
    }
}

template<class Input>
inline void loop_lower_half(int lower_bound, int upper_bound, int inside_diag_num,
                            int left_edge, int top_edge,
                            Input braid_ones,
                            Input mask_prev,
                            Input mask_prev_r,
                            Input *bitset_left_strand_map,
                            Input *bitset_top_strand_map,
                            Input *a_reverse, Input *b) {

    for (int j = lower_bound; j < upper_bound; ++j) {

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
            combing_condition = mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }
    }
}

template<class Input>
inline void loop_upper_half_mpi(int lower_bound, int upper_bound, int rev_counter, int left_edge, int top_edge,
                                Input braid_ones,
                                Input mask,
                                Input mask_r,
                                Input *bitset_left_strand_map,
                                Input *bitset_top_strand_map,
                                Input *a_reverse, Input *b) {

#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
    for (int j = lower_bound; j < upper_bound; ++j) {
        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input left_cap = left_strand >> rev_counter;
        Input symbol_a = a_reverse[left_edge + j];
        Input symbol_b = b[top_edge + j];

        Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        Input combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = combing_condition ^braid_ones;

        if (combing_condition) {
            bitset_top_strand_map[top_edge + j] = (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            top_strand = top_strand << rev_counter;
            symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = mask_r & (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }
    }
}

template<class Input>
inline void loop_center_half_mpi(int lower_bound, int upper_bound,
                                 int left_edge, int top_edge,
                                 Input braid_ones,
                                 Input *bitset_left_strand_map,
                                 Input *bitset_top_strand_map,
                                 Input *a_reverse, Input *b) {
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
    for (int j = lower_bound; j < upper_bound; ++j) {
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

template<class Input>
inline void loop_lower_half_mpi(int lower_bound, int upper_bound, int inside_diag_num,
                                int left_edge, int top_edge,
                                Input braid_ones,
                                Input mask_prev,
                                Input mask_prev_r,
                                Input *bitset_left_strand_map,
                                Input *bitset_top_strand_map,
                                Input *a_reverse, Input *b) {

#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
    for (int j = lower_bound; j < upper_bound; ++j) {

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
            combing_condition = mask_prev_r & (symbols | ((~(left_strand) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }
    }
}


//2 bit per symbol
// a<=b !
template<class Input>
int prefix_lcs_via_braid_bits_4symbol_splited(Input *a_reverse, int a_size, int a_total_symbols,
                                              Input *b, int b_size, int b_total_symbols) {
    int active_symbols_a_active = a_total_symbols % (sizeof(Input) * 8 / 2);
    int active_symbols_b_active = b_total_symbols % (sizeof(Input) * 8 / 2);

    //case 1 x 1
    if (a_size == 1 && b_size == 1 && (active_symbols_a_active != 0 || active_symbols_b_active != 0))
        return prefix_lcs_via_braid_4symbol_one_one_size(a_reverse[0], a_total_symbols, b[0], b_total_symbols);

    // case m x n and one of string not multiple


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



    // phase 0
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {
        loop_upper_half(0, 1, rev_counter, m - 1, 0,
                        braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);
    }

    //PHASE 1
    for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

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


            loop_upper_half(0, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                            braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                            a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half(0, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map, a_reverse, b);

    }

    //PHASE 2
    for (int k = 0; k < total_same_length_diag; ++k) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            loop_upper_half(0, m, rev_counter, 0, k + 1, braid_ones, mask, mask_r, bitset_left_strand_map,
                            bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, m, inside_diag_num, 0, k, braid_ones, mask_prev, mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }
        loop_center_half(0, m, 0, k + 1, braid_ones, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
    }

    //PHASE 3
    auto start_j = total_same_length_diag + 1;
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

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

            loop_upper_half(0, cur_diag_cube_len, rev_counter, 0, start_j, braid_ones, mask, mask_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
            loop_lower_half(0, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, braid_ones, mask_prev,
                            mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half(0, cur_diag_cube_len, 0, start_j, braid_ones, bitset_left_strand_map, bitset_top_strand_map,
                         a_reverse, b);
    }

    //process last triangle in position  m-1, n-1  cube
    mask = braid_ones;
    mask_r = mask;

    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        loop_lower_half(0, 1, inside_diag_num, 0, n - 1, braid_ones, mask, mask_r, bitset_left_strand_map,
                        bitset_top_strand_map,
                        a_reverse, b);
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
int prefix_lcs_via_braid_bits_4symbol_splited_mpi(Input *a_reverse, int a_size, int a_total_symbols,
                                                  Input *b, int b_size, int b_total_symbols, int threads_num) {


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


        Input mask;

#pragma omp  for simd schedule(static)  aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < n; ++k) {
            bitset_top_strand_map[k] = Input(0);
        }

#pragma omp  for simd schedule(static) aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < m; ++k) {
            bitset_left_strand_map[k] = braid_ones;
        }


        auto upper_bound = (sizeof(Input) * 8 / 2) - 1;
        //process first triangle in 0,0 cube
        mask = Input(1);
        Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
        int rev_counter = (sizeof(Input) * 8 - 2);



        // phase 0
#pragma omp single
        {
            for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {
                loop_upper_half(0, 1, rev_counter, m - 1, 0,
                                braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);
            }
        }
        //PHASE 1
        for (int cur_diag_cube_len = 0; cur_diag_cube_len < m; cur_diag_cube_len++) {

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


                loop_upper_half_mpi(0, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                                    braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse,
                                    b);

                loop_lower_half_mpi(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                                    braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                                    a_reverse, b);

                //update mask of current move
                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);
            }

            loop_center_half_mpi(0, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, braid_ones,
                                 bitset_left_strand_map,
                                 bitset_top_strand_map, a_reverse, b);

        }

        //PHASE 2
        for (int k = 0; k < total_same_length_diag; ++k) {

            //to process current
            rev_counter = (sizeof(Input) * 8 - 2);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = braid_ones;
            Input mask_prev_r = braid_ones;


            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
                //update mask of prev move
                mask_prev <<= 2;
                mask_prev_r >>= 2;

                loop_upper_half_mpi(0, m, rev_counter, 0, k + 1, braid_ones, mask, mask_r, bitset_left_strand_map,
                                    bitset_top_strand_map, a_reverse, b);

                loop_lower_half_mpi(0, m, inside_diag_num, 0, k, braid_ones, mask_prev, mask_prev_r,
                                    bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

                //update mask of current move
                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);
            }
            loop_center_half_mpi(0, m, 0, k + 1, braid_ones, bitset_left_strand_map, bitset_top_strand_map, a_reverse,
                                 b);
        }

        //PHASE 3
        auto start_j = total_same_length_diag + 1;
        for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

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

                loop_upper_half_mpi(0, cur_diag_cube_len, rev_counter, 0, start_j, braid_ones, mask, mask_r,
                                    bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
                loop_lower_half_mpi(0, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, braid_ones, mask_prev,
                                    mask_prev_r,
                                    bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

                //update mask of current move
                mask = (mask << 2) | Input(1);
                mask_r = mask_r | (mask_r >> 2);
            }

            loop_center_half_mpi(0, cur_diag_cube_len, 0, start_j, braid_ones, bitset_left_strand_map,
                                 bitset_top_strand_map,
                                 a_reverse, b);
        }

        //process last triangle in position  m-1, n-1  cube
        mask = braid_ones;
        mask_r = mask;

#pragma omp single
        {
            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
                mask = mask << 2;
                mask_r = mask_r >> 2;

                loop_lower_half(0, 1, inside_diag_num, 0, n - 1, braid_ones, mask, mask_r, bitset_left_strand_map,
                                bitset_top_strand_map,
                                a_reverse, b);
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


template<class Input>
int prefix_lcs_via_braid_bits_4symbol_splited_arbitrary(Input *a_reverse, int a_size, int a_total_symbols,
                                                        Input *b, int b_size, int b_total_symbols) {
    int active_symbols_a_active = a_total_symbols % (sizeof(Input) * 8 / 2);
    int active_symbols_b_active = b_total_symbols % (sizeof(Input) * 8 / 2);
//not forget offset!
    //case 1 x 1
    if (a_size == 1 && b_size == 1 && (active_symbols_a_active != 0 || active_symbols_b_active != 0))
        return prefix_lcs_via_braid_4symbol_one_one_size(a_reverse[0], a_total_symbols, b[0], b_total_symbols);


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


    Input l_active_mask = active_symbols_a_active ? Input(0) : braid_ones;
    Input r_active_mask = active_symbols_b_active ? Input(0) : braid_ones;
    for (int i = 0; i < active_symbols_a_active; ++i) {
        l_active_mask |= (Input(1) << ((sizeof(Input) * 8 - 2) - i * 2));
    }
    for (int i = 0; i < active_symbols_b_active; ++i) {
        r_active_mask |= (Input(1) << (2 * i));
    }
    Input tmp_l_active_mask = l_active_mask;
    Input tmp_r_active_mask = r_active_mask;
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


    // handle case when 1xn
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {
        if (a_size == 1) {
            loop_upper_half(0, 1, rev_counter, m - 1, 0, braid_ones, mask & (l_active_mask >> rev_counter),
                            l_active_mask & mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        } else {
            loop_upper_half(0, 1, rev_counter, m - 1, 0, braid_ones, mask, mask_r, bitset_left_strand_map,
                            bitset_top_strand_map, a_reverse, b);
        }
        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);
    }

    //PHASE 1
    for (int cur_diag_cube_len = 0; cur_diag_cube_len < m - 1; cur_diag_cube_len++) {

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


            loop_upper_half(0, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                            braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                            a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half(0, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map, a_reverse, b);

    }


    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

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
//todo
            //current should be splited
            loop_upper_half(0, 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask & (l_active_mask >> rev_counter),
                            mask_r & l_active_mask, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //if squared then  have on both ends
            if (a_size == b_size) {
                loop_upper_half(cur_diag_cube_len, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                                braid_ones, mask & r_active_mask,
                                mask_r & (r_active_mask << rev_counter), bitset_left_strand_map, bitset_top_strand_map,
                                a_reverse, b);
            } else {
                loop_upper_half(cur_diag_cube_len, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                                braid_ones, mask,
                                mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
            }

            loop_upper_half(1, cur_diag_cube_len, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            //previous ok
            loop_lower_half(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                            braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                            a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }


        //first cell
        loop_center_half_with_mask(0, 1, m - 1 - cur_diag_cube_len, 0,
                                   braid_ones, (l_active_mask >> rev_counter), bitset_left_strand_map,
                                   bitset_top_strand_map, a_reverse, b);


        //last cell
        if (a_size == b_size) {
            loop_center_half_with_mask(cur_diag_cube_len, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0,
                                       braid_ones, r_active_mask,
                                       bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        } else {
            loop_center_half_with_mask(cur_diag_cube_len, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0,
                                       braid_ones, braid_ones,
                                       bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        }

        loop_center_half(0, cur_diag_cube_len, m - 1 - cur_diag_cube_len, 0, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map, a_reverse, b);

    }


    //PHASE 2
    for (int k = 0; k < total_same_length_diag; ++k) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            loop_upper_half(0, m, rev_counter, 0, k + 1, braid_ones, mask, mask_r, bitset_left_strand_map,
                            bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, m, inside_diag_num, 0, k, braid_ones, mask_prev, mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }
        loop_center_half(0, m, 0, k + 1, braid_ones, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
    }

    //PHASE 3
    auto start_j = total_same_length_diag + 1;
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

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

            loop_upper_half(0, cur_diag_cube_len, rev_counter, 0, start_j, braid_ones, mask, mask_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
            loop_lower_half(0, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, braid_ones, mask_prev,
                            mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half(0, cur_diag_cube_len, 0, start_j, braid_ones, bitset_left_strand_map, bitset_top_strand_map,
                         a_reverse, b);
    }

    //process last triangle in position  m-1, n-1  cube
    mask = braid_ones;
    mask_r = mask;

    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        loop_lower_half(0, 1, inside_diag_num, 0, n - 1, braid_ones, mask, mask_r, bitset_left_strand_map,
                        bitset_top_strand_map,
                        a_reverse, b);
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
int prefix_lcs_via_braid_bits_4symbol_splited_a_equals_b(Input *a_reverse, int a_size, int a_total_symbols,
                                                         Input *b, int b_size, int b_total_symbols) {

    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;

    Input braid_ones = Input(1);
    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << shift * 2);
    }
    int active_symbols_a_active = a_total_symbols % (sizeof(Input) * 8 / 2);
    int active_symbols_b_active = b_total_symbols % (sizeof(Input) * 8 / 2);

    Input l_active_mask = active_symbols_a_active ? Input(0) : braid_ones;
    Input r_active_mask = active_symbols_b_active ? Input(0) : braid_ones;
    for (int i = 0; i < active_symbols_a_active; ++i) {
        l_active_mask |= (Input(1) << ((sizeof(Input) * 8 - 2) - i * 2));
    }
    for (int i = 0; i < active_symbols_b_active; ++i) {
        r_active_mask |= (Input(1) << (2 * i));
    }

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



    // phase 0
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {
        loop_upper_half(0, 1, rev_counter, m - 1, 0,
                        braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);
    }

    //PHASE 1

    for (int cur_diag_cube_len = 1; cur_diag_cube_len < m - 1; cur_diag_cube_len++) {

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


            loop_upper_half(0, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                            braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                            a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half(0, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map, a_reverse, b);

    }

    //last iter of phase one
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

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


            loop_upper_half(0, 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, Input(mask & (l_active_mask >> rev_counter)), Input(mask_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(cur_diag_cube_len, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, Input(mask & r_active_mask), Input(mask_r & (r_active_mask << rev_counter)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(1, cur_diag_cube_len, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
//tod top edge
            loop_lower_half(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                            braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                            a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half_with_mask(0, 1, m - 1 - cur_diag_cube_len, 0, braid_ones, l_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        loop_center_half_with_mask(cur_diag_cube_len, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, braid_ones,
                                   r_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        loop_center_half(1, cur_diag_cube_len, m - 1 - cur_diag_cube_len, 0, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map, a_reverse, b);
    }



    //PHASE 3
    auto start_j = 0 + 1;
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 2; cur_diag_cube_len--, start_j++) {

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

            loop_upper_half(0, 1, rev_counter, 0, start_j, braid_ones, Input(mask & (l_active_mask >> rev_counter)),
                            Input(mask_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(cur_diag_cube_len - 1, cur_diag_cube_len, rev_counter, 0, start_j, braid_ones,
                            Input(mask & r_active_mask),
                            Input(mask_r & (r_active_mask << rev_counter)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(1, cur_diag_cube_len - 1, rev_counter, 0, start_j, braid_ones, mask, mask_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            loop_lower_half(0, 1, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
            loop_lower_half(cur_diag_cube_len, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & r_active_mask),
                            Input(mask_prev_r & (r_active_mask >> (2 * (inside_diag_num + 1)))),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(1, cur_diag_cube_len, inside_diag_num, 0, start_j - 1, braid_ones, mask_prev, mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }


        loop_center_half_with_mask(0, 1, 0, start_j, braid_ones, l_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        loop_center_half_with_mask(cur_diag_cube_len - 1, cur_diag_cube_len, 0, start_j, braid_ones, r_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        loop_center_half(1, cur_diag_cube_len - 1, 0, start_j, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map,
                         a_reverse, b);
    }


    for (int cur_diag_cube_len = 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

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

            loop_upper_half(0, cur_diag_cube_len, rev_counter, 0, start_j, braid_ones,
                            Input(mask & r_active_mask & (l_active_mask >> rev_counter)),
                            Input(mask_r & l_active_mask & (r_active_mask << rev_counter)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, 1, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & (l_active_mask)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(1, 2, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & r_active_mask),
                            Input(mask_prev_r & (r_active_mask >> (2 * (inside_diag_num + 1)))),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half_with_mask(0, cur_diag_cube_len, 0, start_j, braid_ones, Input(l_active_mask & r_active_mask),
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
    }



    //process last triangle in position  m-1, n-1  cube
    mask = braid_ones;
    mask_r = mask;
    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        loop_lower_half(0, 1, inside_diag_num, 0, n - 1, braid_ones,
                        Input(mask & (l_active_mask << (2 * (inside_diag_num + 1))) & r_active_mask),
                        Input(mask_r & l_active_mask & (r_active_mask >> (2 * (inside_diag_num + 1)))),
                        bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
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


    if(active_symbols_a_active==0) return  a_total_symbols - dis_braid;
    return a_total_symbols - dis_braid + (sizeof(Input) * 4 - active_symbols_a_active);

}


template<class Input>
int prefix_lcs_via_braid_bits_4symbol_splited_a_less_b(Input *a_reverse, int a_size, int a_total_symbols,
                                                       Input *b, int b_size, int b_total_symbols) {

    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;

    Input braid_ones = Input(1);
    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << shift * 2);
    }
    int active_symbols_a_active = a_total_symbols % (sizeof(Input) * 8 / 2);
    int active_symbols_b_active = b_total_symbols % (sizeof(Input) * 8 / 2);

    Input l_active_mask = active_symbols_a_active ? Input(0) : braid_ones;
    Input r_active_mask = active_symbols_b_active ? Input(0) : braid_ones;
    for (int i = 0; i < active_symbols_a_active; ++i) {
        l_active_mask |= (Input(1) << ((sizeof(Input) * 8 - 2) - i * 2));
    }
    for (int i = 0; i < active_symbols_b_active; ++i) {
        r_active_mask |= (Input(1) << (2 * i));
    }

    Input mask;


    for (int k = 0; k < m; ++k) {
        bitset_left_strand_map[k] = braid_ones;
    }

    for (int k = 0; k < n; ++k) {
        bitset_top_strand_map[k] = Input(0);
    }


    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);


    auto upper_bound = (sizeof(Input) * 8 / 2) - 1;
    //process first triangle in 0,0 cube
    mask = Input(1);
    Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
    int rev_counter = (sizeof(Input) * 8 - 2);



    // phase 0
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {
        loop_upper_half(0, 1, rev_counter, m - 1, 0,
                        braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);
    }


    for (int cur_diag_cube_len = 1; cur_diag_cube_len < m - 1; cur_diag_cube_len++) {

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


            loop_upper_half(0, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                            braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                            a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half(0, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map, a_reverse, b);

    }

    //last iter of phase one
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

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

            loop_upper_half(0, 1, rev_counter, 0, 0,
                            braid_ones, Input(mask & (l_active_mask >> rev_counter)), Input(mask_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(1, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0,
                            braid_ones, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                            braid_ones, mask_prev, mask_prev_r, bitset_left_strand_map, bitset_top_strand_map,
                            a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half_with_mask(0, 1, 0, 0, braid_ones, l_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


        loop_center_half(1, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0,
                         braid_ones, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
    }




//todo k-2?вв
    //PHASE 2
    for (int k = 0; k < total_same_length_diag - 1; ++k) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            loop_upper_half(0, 1, rev_counter, 0, k + 1, braid_ones,
                            Input(mask & (l_active_mask >> rev_counter)),
                            Input(mask_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            loop_upper_half(1, m, rev_counter, 0, k + 1, braid_ones,
                            mask, mask_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            loop_lower_half(0, 1, inside_diag_num, 0, k, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & (l_active_mask)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(1, m, inside_diag_num, 0, k, braid_ones,
                            mask_prev,
                            mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half_with_mask(0, 1, 0, k + 1, braid_ones, l_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        loop_center_half(1, m, 0, k + 1, braid_ones, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

    }



    //zero balo
    //last iter of phase 2 both ends
    for (int k = std::max(0, total_same_length_diag - 1); k < total_same_length_diag; ++k) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            loop_upper_half(0, 1, rev_counter, 0, k + 1, braid_ones,
                            Input(mask & (l_active_mask >> rev_counter)),
                            Input(mask_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(m - 1, m, rev_counter, 0, k + 1, braid_ones,
                            Input(mask & r_active_mask),
                            Input(mask_r & (r_active_mask << rev_counter)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            loop_upper_half(1, m - 1, rev_counter, 0, k + 1, braid_ones,
                            mask, mask_r, bitset_left_strand_map,
                            bitset_top_strand_map, a_reverse, b);


            loop_lower_half(0, 1, inside_diag_num, 0, k, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & (l_active_mask)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(1, m, inside_diag_num, 0, k, braid_ones,mask_prev,mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half_with_mask(0, 1, 0, k + 1, braid_ones,
                                   l_active_mask, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        loop_center_half_with_mask(m - 1, m, 0, k + 1, braid_ones,
                                   r_active_mask, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        loop_center_half(1, m - 1, 0, k + 1, braid_ones, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

    }


    //phase 3
    auto start_j = total_same_length_diag + 1;
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 2; cur_diag_cube_len--, start_j++) {

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

            loop_upper_half(0, 1, rev_counter, 0, start_j, braid_ones, Input(mask & (l_active_mask >> rev_counter)),
                            Input(mask_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(cur_diag_cube_len - 1, cur_diag_cube_len, rev_counter, 0, start_j, braid_ones,
                            Input(mask & r_active_mask),
                            Input(mask_r & (r_active_mask << rev_counter)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_upper_half(1, cur_diag_cube_len - 1, rev_counter, 0, start_j, braid_ones, mask, mask_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            loop_lower_half(0, 1, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
            loop_lower_half(cur_diag_cube_len, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & r_active_mask),
                            Input(mask_prev_r & (r_active_mask >> (2 * (inside_diag_num + 1)))),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(1, cur_diag_cube_len, inside_diag_num, 0, start_j - 1, braid_ones, mask_prev, mask_prev_r,
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }


        loop_center_half_with_mask(0, 1, 0, start_j, braid_ones, l_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        loop_center_half_with_mask(cur_diag_cube_len - 1, cur_diag_cube_len, 0, start_j, braid_ones, r_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        loop_center_half(1, cur_diag_cube_len - 1, 0, start_j, braid_ones, bitset_left_strand_map,
                         bitset_top_strand_map, a_reverse, b);
    }




//todo
    for (int cur_diag_cube_len = 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

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

            loop_upper_half(0, 1, rev_counter, 0, start_j, braid_ones,
                            Input(mask & r_active_mask & (l_active_mask >> rev_counter)),
                            Input(mask_r & l_active_mask & (r_active_mask << rev_counter)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(0, 1, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & (l_active_mask)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half(1, 2, inside_diag_num, 0, start_j - 1, braid_ones,
                            Input(mask_prev & r_active_mask),
                            Input(mask_prev_r & (r_active_mask >> (2 * (inside_diag_num + 1)))),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        // todo no need?
        loop_center_half_with_mask(0, 1, 0, start_j, braid_ones, Input(l_active_mask & r_active_mask),
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
    }



    //last
    //process last triangle in position  m-1, n-1  cube
    mask = braid_ones;
    mask_r = mask;
    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        loop_lower_half(0, 1, inside_diag_num, 0, n - 1, braid_ones,
                        Input(mask & (l_active_mask << (2 * (inside_diag_num + 1))) & r_active_mask),
                        Input(mask_r & l_active_mask & (r_active_mask >> (2 * (inside_diag_num + 1)))),
                        bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
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


    if(active_symbols_a_active==0) return  a_total_symbols - dis_braid;
    return a_total_symbols - dis_braid + (sizeof(Input) * 4 - active_symbols_a_active);

}


template<class Input>
int prefix_lcs_via_braid_bits_4symbol_splited_a_one_b_arbitrary(Input *a_reverse, int a_size, int a_total_symbols,
                                                       Input *b, int b_size, int b_total_symbols) {

    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;

    Input braid_ones = Input(1);
    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << shift * 2);
    }


    int active_symbols_a_active = a_total_symbols % (sizeof(Input) * 8 / 2);
    int active_symbols_b_active = b_total_symbols % (sizeof(Input) * 8 / 2);

    Input l_active_mask = active_symbols_a_active ? Input(0) : braid_ones;
    Input r_active_mask = active_symbols_b_active ? Input(0) : braid_ones;
    for (int i = 0; i < active_symbols_a_active; ++i) {
        l_active_mask |= (Input(1) << ((sizeof(Input) * 8 - 2) - i * 2));
    }
    for (int i = 0; i < active_symbols_b_active; ++i) {
        r_active_mask |= (Input(1) << (2 * i));
    }

    Input mask;


    for (int k = 0; k < m; ++k) {
        bitset_left_strand_map[k] = braid_ones;
    }

    for (int k = 0; k < n; ++k) {
        bitset_top_strand_map[k] = Input(0);
    }


    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);


    auto upper_bound = (sizeof(Input) * 8 / 2) - 1;
    //process first triangle in 0,0 cube
    mask = Input(1);
    Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
    int rev_counter = (sizeof(Input) * 8 - 2);



    // phase 0
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter -= 2) {
        loop_upper_half(0, 1, rev_counter, m - 1, 0,
                        braid_ones, Input(mask & (l_active_mask>>rev_counter)), Input(mask_r& l_active_mask),
                        bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);
    }



    //PHASE 2
    for (int k = 0; k < total_same_length_diag - 1; ++k) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            loop_upper_half(0, 1, rev_counter, 0, k + 1, braid_ones,
                            Input(mask & (l_active_mask >> rev_counter)),
                            Input(mask_r & l_active_mask),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            loop_lower_half(0, 1, inside_diag_num, 0, k, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & (l_active_mask)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half_with_mask(0, 1, 0, k + 1, braid_ones, l_active_mask,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

    }

    for (int k = std::max(0, total_same_length_diag - 1); k < total_same_length_diag; ++k) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 2);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = braid_ones;
        Input mask_prev_r = braid_ones;


        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
            //update mask of prev move
            mask_prev <<= 2;
            mask_prev_r >>= 2;

            loop_upper_half(0, 1, rev_counter, 0, k + 1, braid_ones,
                            Input(mask & (l_active_mask >> rev_counter) & r_active_mask),
                            Input(mask_r & l_active_mask & (r_active_mask<<rev_counter)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);


            loop_lower_half(0, 1, inside_diag_num, 0, k, braid_ones,
                            Input(mask_prev & (l_active_mask << (2 * (inside_diag_num + 1)))),
                            Input(mask_prev_r & (l_active_mask)),
                            bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 2) | Input(1);
            mask_r = mask_r | (mask_r >> 2);
        }

        loop_center_half_with_mask(0, 1, 0, k + 1, braid_ones,
                                   l_active_mask & r_active_mask, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

    }


    //last
    //process last triangle in position  m-1, n-1  cube
    mask = braid_ones;
    mask_r = mask;
    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 2;
        mask_r = mask_r >> 2;

        loop_lower_half(0, 1, inside_diag_num, 0, n - 1, braid_ones,
                        Input(mask & (l_active_mask << (2 * (inside_diag_num + 1))) & r_active_mask),
                        Input(mask_r & l_active_mask & (r_active_mask >> (2 * (inside_diag_num + 1)))),
                        bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
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


    if(active_symbols_a_active==0) return  a_total_symbols - dis_braid;
    return a_total_symbols - dis_braid + (sizeof(Input) * 4 - active_symbols_a_active);

}




#endif //CPU_TMP_H
