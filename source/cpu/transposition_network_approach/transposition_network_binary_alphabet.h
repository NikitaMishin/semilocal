//
// Created by nikita on 19.08.2020.
//

#ifndef CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H
#define CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H


#include <stdlib.h>


template<class Input>
inline void loop_upper_half_binary(int lower_bound, int upper_bound, int rev_counter, int left_edge, int top_edge,
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
        Input combing_condition =
                mask & (((~(symbol_a >> rev_counter)) ^ symbol_b) | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = ~combing_condition;

        if (combing_condition) {
            bitset_top_strand_map[top_edge + j] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            top_strand = top_strand << rev_counter;
            combing_condition =
                    mask_r &
                    (((~(symbol_a)) ^ (symbol_b << rev_counter)) | ((~(left_strand) & top_strand)));
            rev_combing_cond = ~combing_condition;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }
    }
}

template<class Input>
inline void loop_center_half_binary(int lower_bound, int upper_bound,
                                    int left_edge, int top_edge,
                                    Input *bitset_left_strand_map,
                                    Input *bitset_top_strand_map,
                                    Input *a_reverse, Input *b) {

    for (int j = lower_bound; j < upper_bound; ++j) {
        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) | ((~left_strand) & top_strand));
        Input rev_combing_cond = ~combing_condition;
        if (combing_condition) {
            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
            bitset_top_strand_map[top_edge + j] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_strand);
        }
    }
}


template<class Input>
inline void loop_lower_half_binary(int lower_bound, int upper_bound, int inside_diag_num,
                                   int left_edge, int top_edge,
                                   Input mask_prev,
                                   Input mask_prev_r,
                                   Input *bitset_left_strand_map,
                                   Input *bitset_top_strand_map,
                                   Input *a_reverse, Input *b) {

    for (int j = lower_bound; j < upper_bound; ++j) {

        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input left_cap = left_strand << (inside_diag_num + 1);
        Input symbol_a = a_reverse[left_edge + j];
        Input symbol_b = b[top_edge + j];
        Input combing_condition =
                mask_prev &
                (((~(symbol_a << (inside_diag_num + 1))) ^ symbol_b) | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = ~combing_condition;

        if (combing_condition) {
            bitset_top_strand_map[top_edge + j] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            top_strand = top_strand >> (inside_diag_num + 1);
            combing_condition =
                    mask_prev_r &
                    (((~(symbol_a)) ^ (symbol_b >> (inside_diag_num + 1))) |
                     ((~(left_strand) & top_strand)));
            rev_combing_cond = ~combing_condition;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }

    }
}

template<class Input>
inline void loop_upper_half_mpi_binary(int lower_bound, int upper_bound, int rev_counter, int left_edge, int top_edge,
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
        Input combing_condition =
                mask & (((~(symbol_a >> rev_counter)) ^ symbol_b) | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = ~combing_condition;

        if (combing_condition) {
            bitset_top_strand_map[top_edge + j] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            top_strand = top_strand << rev_counter;
            combing_condition =
                    mask_r &
                    (((~(symbol_a)) ^ (symbol_b << rev_counter)) | ((~(left_strand) & top_strand)));
            rev_combing_cond = ~combing_condition;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }
    }
}

template<class Input>
inline void loop_center_half_mpi_binary(int lower_bound, int upper_bound,
                                        int left_edge, int top_edge,
                                        Input *bitset_left_strand_map,
                                        Input *bitset_top_strand_map,
                                        Input *a_reverse, Input *b) {

#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
    for (int j = lower_bound; j < upper_bound; ++j) {
        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) | ((~left_strand) & top_strand));
        Input rev_combing_cond = ~combing_condition;
        if (combing_condition) {
            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
            bitset_top_strand_map[top_edge + j] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_strand);
        }
    }
}

template<class Input>
inline void loop_lower_half_mpi_binary(int lower_bound, int upper_bound, int inside_diag_num,
                                       int left_edge, int top_edge,
                                       Input mask_prev,
                                       Input mask_prev_r,
                                       Input *bitset_left_strand_map,
                                       Input *bitset_top_strand_map,
                                       Input *a_reverse, Input *b) {

#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
    for (int j = lower_bound; j < upper_bound; ++j) {

        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];
        Input left_cap = left_strand << (inside_diag_num + 1);
        Input symbol_a = a_reverse[left_edge + j];
        Input symbol_b = b[top_edge + j];
        Input combing_condition =
                mask_prev &
                (((~(symbol_a << (inside_diag_num + 1))) ^ symbol_b) | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = ~combing_condition;

        if (combing_condition) {
            bitset_top_strand_map[top_edge + j] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);
            top_strand = top_strand >> (inside_diag_num + 1);
            combing_condition =
                    mask_prev_r &
                    (((~(symbol_a)) ^ (symbol_b >> (inside_diag_num + 1))) |
                     ((~(left_strand) & top_strand)));
            rev_combing_cond = ~combing_condition;

            bitset_left_strand_map[left_edge + j] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }
    }
}


////assume both multiple by sizeof(Input)*8
template<class Input>
int prefix_lcs_via_braid_bits_binary(Input *a_reverse, int a_size, int a_total_symbols,
                                     Input *b, int b_size, int b_total_symbols) {


    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);


    Input mask;


    for (int k = 0; k < m; ++k) {
        bitset_left_strand_map[k] = ~Input(0);
    }

    for (int k = 0; k < n; ++k) {
        bitset_top_strand_map[k] = Input(0);
    }

    auto upper_bound = (sizeof(Input) * 8) - 1;
    //process first triangle in 0,0 cube
    mask = Input(1);
    Input mask_r = Input(1) << (sizeof(Input) * 8 - 1);
    int rev_counter = (sizeof(Input) * 8 - 1);

    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter--) {
        loop_upper_half_binary(0, 1, rev_counter, m - 1, 0, mask, mask_r, bitset_left_strand_map, bitset_top_strand_map,
                               a_reverse, b);
        mask = (mask << 1) | Input(1);
        mask_r = mask_r | (mask_r >> 1);
    }

    //PHASE 1: Process diagonal till fill big left triangle,
    for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {
        //to process current
        rev_counter = (sizeof(Input) * 8 - 1);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = ~static_cast<Input>(0);
        Input mask_prev_r = ~Input(0);


        //process previous size/2 - 1 cubes and current size/2 -1  cubes
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
            //update mask of prev move
            mask_prev <<= 1;
            mask_prev_r >>= 1;


            loop_upper_half_binary(0, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0, mask, mask_r,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            loop_lower_half_binary(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0, mask_prev,
                                   mask_prev_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
            //update mask of current move
            mask = (mask << 1) | Input(1);
            mask_r = mask_r | (mask_r >> 1);
        }

        loop_center_half_binary(0, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, bitset_left_strand_map,
                                bitset_top_strand_map, a_reverse, b);
    }

    //PHASE 2
    for (int k = 0; k < total_same_length_diag; ++k) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 1);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = ~Input(0);
        Input mask_prev_r = ~Input(0);

        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
            //update mask of prev move
            mask_prev <<= 1;
            mask_prev_r >>= 1;


            loop_upper_half_binary(0, m, rev_counter, 0, k + 1, mask, mask_r, bitset_left_strand_map,
                                   bitset_top_strand_map,
                                   a_reverse, b);

            loop_lower_half_binary(0, m, inside_diag_num, 0, k, mask_prev, mask_prev_r, bitset_left_strand_map,
                                   bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 1) | Input(1);
            mask_r = mask_r | (mask_r >> 1);
        }

        loop_center_half_binary(0, m, 0, k + 1, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
    }


    auto start_j = total_same_length_diag + 1;
    for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

        //to process current
        rev_counter = (sizeof(Input) * 8 - 1);
        mask = Input(1);
        mask_r = Input(1) << rev_counter;

        //to process previous
        Input mask_prev = ~Input(0);
        Input mask_prev_r = ~Input(0);

        //process previous size/2 - 1 cubes and current size/2 -1  cubes
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
            //update mask of prev move
            mask_prev <<= 1;
            mask_prev_r >>= 1;

            loop_upper_half_binary(0, cur_diag_cube_len, rev_counter, 0, start_j, mask, mask_r, bitset_left_strand_map,
                                   bitset_top_strand_map, a_reverse, b);

            loop_lower_half_binary(0, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, mask_prev, mask_prev_r,
                                   bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

            //update mask of current move
            mask = (mask << 1) | Input(1);
            mask_r = mask_r | (mask_r >> 1);

        }

        loop_center_half_binary(0, cur_diag_cube_len, 0, start_j, bitset_left_strand_map, bitset_top_strand_map,
                                a_reverse, b);
    }


    //process last triangle in position  m-1, n-1  cube
    mask = ~Input(0);
    mask_r = mask;

    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 1;
        mask_r = mask_r >> 1;


        loop_lower_half_binary(0, 1, inside_diag_num, 0, n - 1, mask, mask_r, bitset_left_strand_map,
                               bitset_top_strand_map, a_reverse, b);
    }

    for (int i1 = 0; i1 < m; ++i1) {
        //  Brian Kernighan’s Algorithm
        auto counter = 0;
        auto number = bitset_left_strand_map[i1];
        //  LogNumber
        while (number) {
            number &= (number - 1);
            counter++;
        }
        dis_braid += counter;
    }

    free(bitset_top_strand_map);
    free(bitset_left_strand_map);


    return a_total_symbols - dis_braid;

}


//assume both multiple by sizeof(Input)*8
template<class Input>
int prefix_lcs_via_braid_bits_binary_mpi(Input *a_reverse, int a_size, int a_total_symbols,
                                         Input *b, int b_size, int b_total_symbols,
                                         int threads_num = 1) {


    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);

#pragma omp parallel num_threads(threads_num)  default(none) shared(bitset_left_strand_map, bitset_top_strand_map, a_reverse, b, m, n, dis_braid, total_same_length_diag)
    {

        int left_edge, top_edge;
        Input mask;


#pragma omp  for simd schedule(static) aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < m; ++k) {
            bitset_left_strand_map[k] = ~Input(0);
        }

#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map:sizeof(Input)*8)
        for (int k = 0; k < n; ++k) {
            bitset_top_strand_map[k] = Input(0);
        }

        auto upper_bound = (sizeof(Input) * 8) - 1;
        //process first triangle in 0,0 cube
        mask = Input(1);
        Input mask_r = Input(1) << (sizeof(Input) * 8 - 1);
        int rev_counter = (sizeof(Input) * 8 - 1);

        //PHASE 0:Process first triangle
#pragma omp single
        {
            for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter--) {
                loop_upper_half_binary(0, 1, rev_counter, m - 1, 0, mask, mask_r, bitset_left_strand_map,
                                       bitset_top_strand_map, a_reverse, b);
                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);
            }
        }

        //PHASE 1: Process diagonal till fill big left triangle,
        for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {
            //to process current
            rev_counter = (sizeof(Input) * 8 - 1);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = ~static_cast<Input>(0);
            Input mask_prev_r = ~Input(0);


            //process previous size/2 - 1 cubes and current size/2 -1  cubes
            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
                //update mask of prev move
                mask_prev <<= 1;
                mask_prev_r >>= 1;


                loop_upper_half_mpi_binary(0, cur_diag_cube_len + 1, rev_counter, m - 1 - cur_diag_cube_len, 0, mask,
                                           mask_r,
                                           bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

                loop_lower_half_mpi_binary(0, cur_diag_cube_len, inside_diag_num, m - 1 - cur_diag_cube_len + 1, 0,
                                           mask_prev,
                                           mask_prev_r, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
                //update mask of current move
                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);
            }

            loop_center_half_mpi_binary(0, cur_diag_cube_len + 1, m - 1 - cur_diag_cube_len, 0, bitset_left_strand_map,
                                        bitset_top_strand_map, a_reverse, b);
        }

        //PHASE 2
        for (int k = 0; k < total_same_length_diag; ++k) {

            //to process current
            rev_counter = (sizeof(Input) * 8 - 1);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = ~Input(0);
            Input mask_prev_r = ~Input(0);

            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
                //update mask of prev move
                mask_prev <<= 1;
                mask_prev_r >>= 1;


                loop_upper_half_mpi_binary(0, m, rev_counter, 0, k + 1, mask, mask_r, bitset_left_strand_map,
                                           bitset_top_strand_map,
                                           a_reverse, b);

                loop_lower_half_mpi_binary(0, m, inside_diag_num, 0, k, mask_prev, mask_prev_r, bitset_left_strand_map,
                                           bitset_top_strand_map, a_reverse, b);

                //update mask of current move
                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);
            }

            loop_center_half_mpi_binary(0, m, 0, k + 1, bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);
        }


        auto start_j = total_same_length_diag + 1;
        for (int cur_diag_cube_len = m - 1; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {

            //to process current
            rev_counter = (sizeof(Input) * 8 - 1);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = ~Input(0);
            Input mask_prev_r = ~Input(0);

            //process previous size/2 - 1 cubes and current size/2 -1  cubes
            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter--) {
                //update mask of prev move
                mask_prev <<= 1;
                mask_prev_r >>= 1;

                loop_upper_half_mpi_binary(0, cur_diag_cube_len, rev_counter, 0, start_j, mask, mask_r,
                                           bitset_left_strand_map,
                                           bitset_top_strand_map, a_reverse, b);

                loop_lower_half_mpi_binary(0, cur_diag_cube_len + 1, inside_diag_num, 0, start_j - 1, mask_prev,
                                           mask_prev_r,
                                           bitset_left_strand_map, bitset_top_strand_map, a_reverse, b);

                //update mask of current move
                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);

            }

            loop_center_half_mpi_binary(0, cur_diag_cube_len, 0, start_j, bitset_left_strand_map, bitset_top_strand_map,
                                        a_reverse, b);
        }


        //process last triangle in position  m-1, n-1  cube
        mask = ~Input(0);
        mask_r = mask;
#pragma omp single
        {

            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
                mask = mask << 1;
                mask_r = mask_r >> 1;


                loop_lower_half_binary(0, 1, inside_diag_num, 0, n - 1, mask, mask_r, bitset_left_strand_map,
                                       bitset_top_strand_map, a_reverse, b);
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

    return a_total_symbols - dis_braid;

}


#endif //CPU_TRANSPOSITION_NETWORK_BINARY_ALPHABET_H
