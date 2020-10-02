//
// Created by nikita on 30.09.2020.
//

#ifndef CPU_2SYMBOL_NEW_H
#define CPU_2SYMBOL_NEW_H

template<class Input>
void precompute_info_table(int *a, int a_size, int *b, int b_size) {
//    int bit_per_word = 8 * sizeof(Input);
//    for (int i = 0; i < ; ++i) {
//        // load 64 elems
//        // 64 equals
//
//    }
}


#include <vector>
#include <cmath>
#include <bitset>


template<class Input>
inline void process_cubes_antidiag_bin(int lower_bound, int upper_bound, int left_edge, int top_edge,
                                       Input braid_ones,
                                       Input *bitset_left_strand_map,
                                       Input *bitset_top_strand_map,
                                       Input *a_reverse, Input *b) {

    for (int j = lower_bound; j < upper_bound; ++j) {
        Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;

        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input symbol_a = a_reverse[left_edge + j];
        Input symbol_b = b[top_edge + j];

        int rev_counter = (sizeof(Input) * 8 - 2);
        Input mask = Input(1);
//        Input mask_r = Input(1) << rev_counter;


        // upper half
#pragma GCC unroll  128
        for (int inside_diag_num = 0; inside_diag_num < sizeof(Input) * 8 / 2 - 1; ++inside_diag_num) {
            left_cap = left_strand >> rev_counter;
            symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            if (combing_condition) {
                top_strand_shifted = top_strand << rev_counter;
                top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);


//                symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
//                symbols &= (symbols >> 1) & braid_ones;
//                combing_condition = mask_r & (symbols | ((~(left_strand) & top_strand_shifted)));

                combing_condition <<= rev_counter;
                rev_combing_cond = combing_condition ^ braid_ones;

                left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
            }


            rev_counter -= 2;
            mask = (mask << 2) | Input(1);
//            mask_r = mask_r | (mask_r >> 2);
        }

        // center
        symbols = (~(symbol_a ^ symbol_b));
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = (symbols | ((~left_strand) & top_strand));
        rev_combing_cond = combing_condition ^ braid_ones;
        if (combing_condition) {
            top_strand_shifted = top_strand;
            top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
        }


        mask = braid_ones;
//        mask_r = braid_ones;

        //lower half
#pragma GCC unroll  128
        for (int inside_diag_num = 0; inside_diag_num < sizeof(Input) * 8 / 2 - 1; ++inside_diag_num) {
            mask <<= 2;
//            mask_r >>= 2;

            left_cap = left_strand << (2 * (inside_diag_num + 1));
            symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
            symbols &= (symbols >> 1) & braid_ones;

            combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));
            rev_combing_cond = combing_condition ^ braid_ones;

            if (combing_condition) {
                top_strand_shifted = top_strand >> (2 * (inside_diag_num + 1));
                top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);
//                symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
//                symbols &= (symbols >> 1) & braid_ones;

//                combing_condition = mask_r & (symbols | ((~(left_strand) & top_strand_shifted)));
                combing_condition >>= (2 * (inside_diag_num + 1));
                rev_combing_cond = combing_condition ^ braid_ones;


                left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
            }
        }


        bitset_left_strand_map[left_edge + j] = left_strand;

        bitset_top_strand_map[top_edge + j] = top_strand;
    }
}

template<class Input>
inline void process_cubes_antidiag_mpi_bin(int lower_bound, int upper_bound, int left_edge, int top_edge,
                                           Input *bitset_left_strand_map,
                                           Input *bitset_top_strand_map,
                                           Input *a_reverse, Input *b) {

    const int upper = sizeof(Input) * 8 - 1;

#pragma omp   for  simd schedule(static)  aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
    for (int j = lower_bound; j < upper_bound; ++j) {
        Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;
        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input symbol_a = a_reverse[left_edge + j];
        Input symbol_b = b[top_edge + j];

        Input mask = Input(1);


        // upper half
#pragma GCC unroll  256
        for (int rev_counter = (sizeof(Input) * 8 - 1); rev_counter > 0; rev_counter--) {
            left_cap = left_strand >> rev_counter;
            symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
            combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));
            rev_combing_cond = ~combing_condition;

            top_strand_shifted = top_strand << rev_counter;
            top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            combing_condition <<= rev_counter;
            rev_combing_cond = ~combing_condition;

            left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);

            mask = (mask << 1) | Input(1);
        }

        // center
        symbols = (~(symbol_a ^ symbol_b));
//        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = (symbols | ((~left_strand) & top_strand));
        rev_combing_cond = ~combing_condition;
        top_strand_shifted = top_strand;
        top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_strand);
        left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);

        mask = ~Input(0);

        //lower half
#pragma GCC unroll 256
        for (int inside_diag_num = 1; inside_diag_num < upper + 1; inside_diag_num++) {
            mask <<= 1;

            left_cap = left_strand << (inside_diag_num);
            symbols = ~(((symbol_a << inside_diag_num)) ^ symbol_b);
//            symbols &= (symbols >> 1) & braid_ones;

            combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));
            rev_combing_cond = ~combing_condition;

            top_strand_shifted = top_strand >> ((inside_diag_num));
            top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            combing_condition >>= ((inside_diag_num));
            rev_combing_cond = ~combing_condition;


            left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
        }


        bitset_left_strand_map[left_edge + j] = left_strand;
        bitset_top_strand_map[top_edge + j] = top_strand;
    }
}


template<class Input>
inline void process_cube_with_exception_bin(int left_edge, int top_edge, int j, Input braid_ones, Input l_active_mask,
                                            Input r_active_mask,
                                            Input *bitset_left_strand_map, Input *bitset_top_strand_map,
                                            Input *a_reverse,
                                            Input *b) {

    Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;

    Input left_strand = bitset_left_strand_map[left_edge + j];
    Input top_strand = bitset_top_strand_map[top_edge + j];
    Input symbol_a = a_reverse[left_edge + j];
    Input symbol_b = b[top_edge + j];

    int rev_counter = (sizeof(Input) * 8 - 2);
    Input mask = Input(1);
    Input mask_r = Input(1) << rev_counter;


    // upper half
    for (int inside_diag_num = 0; inside_diag_num < sizeof(Input) * 8 / 2; ++inside_diag_num) {
        left_cap = left_strand >> rev_counter;
        symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition =
                r_active_mask & (l_active_mask >> rev_counter) & mask & (symbols | (((~(left_cap)) & top_strand)));
        rev_combing_cond = combing_condition ^ braid_ones;

        if (combing_condition) {
            top_strand_shifted = top_strand << rev_counter;
            top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            symbols = ~(((symbol_a)) ^ (symbol_b << rev_counter));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = (r_active_mask << rev_counter) & l_active_mask & mask_r &
                                (symbols | ((~(left_strand) & top_strand_shifted)));
            rev_combing_cond = combing_condition ^ braid_ones;

            left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
        }


        rev_counter -= 2;
        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);
    }

    // center
    symbols = (~(symbol_a ^ symbol_b));
    symbols &= (symbols >> 1) & braid_ones;
    combing_condition = l_active_mask & r_active_mask & (symbols | ((~left_strand) & top_strand));
    rev_combing_cond = combing_condition ^ braid_ones;
    if (combing_condition) {
        top_strand_shifted = top_strand;
        top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_strand);
        left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
    }


    mask = braid_ones;
    mask_r = braid_ones;

    //lower half
    for (int inside_diag_num = 0; inside_diag_num < sizeof(Input) * 8 / 2; ++inside_diag_num) {
        mask <<= 2;
        mask_r >>= 2;

        left_cap = left_strand << (2 * (inside_diag_num + 1));
        symbols = ~(((symbol_a << (2 * (inside_diag_num + 1)))) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        combing_condition = r_active_mask & (l_active_mask << (2 * (inside_diag_num + 1))) & mask &
                            (symbols | (((~(left_cap)) & top_strand)));
        rev_combing_cond = combing_condition ^ braid_ones;

        if (combing_condition) {
            top_strand_shifted = top_strand >> (2 * (inside_diag_num + 1));
            top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            symbols = ~(((symbol_a)) ^ (symbol_b >> (2 * (inside_diag_num + 1))));
            symbols &= (symbols >> 1) & braid_ones;

            combing_condition = (r_active_mask >> (2 * (inside_diag_num + 1))) & l_active_mask & mask_r &
                                (symbols | ((~(left_strand) & top_strand_shifted)));
            rev_combing_cond = combing_condition ^ braid_ones;

            left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
        }
    }


    bitset_left_strand_map[left_edge + j] = left_strand;
    bitset_top_strand_map[top_edge + j] = top_strand;

}


template<class Input>
int prefix_lcs_via_braid_bits_2symbol_v2_full_mask(Input *a_reverse, int a_size, int a_total_symbols,
                                                   Input *b, int b_size, int b_total_symbols, int threads_num) {


    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));


    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);

    Input braid_ones = ~Input(0);


#pragma omp parallel num_threads(threads_num) default(none) shared(bitset_left_strand_map, bitset_top_strand_map, a_reverse, b, m, n, dis_braid, total_same_length_diag, braid_ones)
    {

#pragma omp  for simd schedule(static)  aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < n; ++k) {
            bitset_top_strand_map[k] = Input(0);
        }

#pragma omp  for simd schedule(static) aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < m; ++k) {
            bitset_left_strand_map[k] = braid_ones;
        }

        for (int diag_len = 0; diag_len < m - 1; diag_len++) {
            process_cubes_antidiag_mpi_bin(0, diag_len + 1, m - 1 - diag_len, 0, bitset_left_strand_map,
                                           bitset_top_strand_map, a_reverse, b);

        }

        for (int k = 0; k < total_same_length_diag; k++) {
            process_cubes_antidiag_mpi_bin(0, m, 0, k, bitset_left_strand_map,
                                           bitset_top_strand_map, a_reverse, b);
        }

        auto start_j = total_same_length_diag;

        for (int diag_len = m - 1; diag_len >= 1; diag_len--) {
            process_cubes_antidiag_mpi_bin(0, diag_len, 0, start_j, bitset_left_strand_map,
                                           bitset_top_strand_map, a_reverse, b);
            start_j++;
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


#endif //CPU_2SYMBOL_NEW_H
