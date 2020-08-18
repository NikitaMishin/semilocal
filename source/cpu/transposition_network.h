//
// Created by nikita on 30.07.2020.
//

#ifndef CPU_TRANSPOSITION_NETWORK_H
#define CPU_TRANSPOSITION_NETWORK_H

#include <vector>
#include <cmath>

template<class Input>
int prefix_lcs_via_braid_sequential(std::vector<Input> const &a, std::vector<Input> const &b) {

    auto m = a.size();
    auto n = b.size();

    auto dis_braid = 0;
    auto size = m + n;
    auto strand_map = new bool[size];

    for (int k = 0; k < m; ++k) {
        strand_map[k] = true;
    }

    for (int l = m; l < m + n; ++l) {
        strand_map[l] = false;
    }


    for (int i = 0; i < m; ++i) {
        for (int j = 0, top_edge = m; j < n; ++j, top_edge++) {
            bool left_strand = strand_map[i];
            bool right_strand = strand_map[top_edge];
            bool r = a[i] == b[j] || (!left_strand && right_strand);
            if (r) std::swap(strand_map[top_edge], strand_map[i]);
//            strand_map[left_edge] = !(!right_strand && (d == left_strand || ((d && !left_strand))));
//            strand_map[left_edge] = !(d || right_strand);
//            strand_map[top_edge] = (d && left_strand) || (left_strand && right_strand);
//            strand_map[top_edge] = (!left_strand && right_strand)? left_strand :    left_strand || right_strand &&d ;
//            strand_map[left_edge] =  right_strand || left_strand || !d ;
        }
    }


    for (int i1 = 0; i1 < m; ++i1) {
        dis_braid += strand_map[i1];
    }


    return m - dis_braid;
}

/**
 * Assume a <= b
 * unsigned types only!
 * @tparam Input
 * @param a
 * @param b
 * @return
 */
template<class Input>
int prefix_lcs_via_braid_mpi(std::vector<Input> const &a, std::vector<Input> const &b, int threads_num = 1) {

    auto m = a.size();
    auto n = b.size();

    auto dis_braid = 0;
    auto size = m + n;
//    short
    auto strand_map = new int[size];

    auto num_diag = a.size() + b.size() - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);


#pragma omp parallel num_threads(threads_num)  default(none) shared(a, b, strand_map, total_same_length_diag, size, m, n, dis_braid)
    {
        int left_edge, top_edge;
        //    init phase
#pragma omp  for simd schedule(static)
        for (int k = 0; k < m; ++k) {
            strand_map[k] = true;
        }

#pragma omp for simd schedule(static)
        for (int l = 0; l < n; ++l) {
            strand_map[l + m] = false;
        }

        //    phase one
        for (int cur_diag_len = 0; cur_diag_len < m - 1; ++cur_diag_len) {
            left_edge = m - 1 - cur_diag_len;
            top_edge = m;
#pragma omp for simd schedule(static)
            for (int j = 0; j < cur_diag_len + 1; ++j) {
                auto left_strand = strand_map[left_edge + j];
                auto right_strand = strand_map[top_edge + j];
                auto r = a[cur_diag_len - j] == b[j] || (!left_strand && right_strand);
//                bool r = a[left_edge + j] == b[j] || (!left_strand && right_strand);
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
//  ||  ||  ||   || auto r = a[left_edge + k] == b[left_edge + j + k] || (!left_strand && right_strand);
//  ||  ||  ||   || auto r = |  |  a[left_edge + k] == b[left_edge + j + k]  |  |
//  ||  ||  ||   ||
                bool r = a[i - k] == b[left_edge + j + k] || (!left_strand && right_strand);
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
                auto left_strand = strand_map[left_edge + k];
                auto right_strand = strand_map[top_edge + k];
//                auto r = a[left_edge+k] == b[j + k] || (!left_strand && right_strand);
                bool r = a[i - k] == b[j + k] || (!left_strand && right_strand);
                if (r) std::swap(strand_map[top_edge + k], strand_map[left_edge + k]);
            }
        }

#pragma omp for simd reduction(+:dis_braid) schedule(static)
        for (int i1 = 0; i1 < m; ++i1) {
            dis_braid += strand_map[i1];
        }

    }

    delete[] strand_map;


    return m - dis_braid;
}


//assume both multiple by sizeof(Input)*8
template<class Input>
int prefix_lcs_via_braid_bits_binary_mpi(Input *a_reverse, int a_size, int a_total_symbols,
                                         Input *b, int b_size, int b_total_symbols, int alphabet_size,
                                         int threads_num = 1) {


    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);

//    auto not_garbage_bits_a = a_total_symbols % (sizeof(Input) * 8);
//    auto not_garbage_bits_b = b_total_symbols % (sizeof(Input) * 8);


//    Input not_garbage_mask_a_rev, not_garbage_mask_b;
//    //all ok
//    if(not_garbage_bits_a==0){
//        not_garbage_mask_a_rev = ~Input(0);
//    } else {
//        not_garbage_mask_a_rev = (~Input(0)) ^ ((Input(1) << (sizeof(Input)*8 - not_garbage_bits_a)) - 1);
//    }
//    if(not_garbage_bits_b==0){
//        not_garbage_mask_b = ~Input(0);
//    } else {
//        not_garbage_mask_b = (Input(1) << not_garbage_bits_b) - 1;
//    }
//,not_garbage_mask_a_rev, not_garbage_mask_b


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

                Input left_strand = bitset_left_strand_map[m - 1];
                Input top_strand = bitset_top_strand_map[0];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[m - 1];
                Input symbol_b = b[0];

                Input combing_condition =
                        mask & (((~(symbol_a >> rev_counter)) ^ symbol_b) | (((~(left_cap)) & top_strand)));
                Input rev_combing_cond = ~combing_condition;

                if (combing_condition) {
                    bitset_top_strand_map[0] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);

                    top_strand = top_strand << rev_counter;
                    combing_condition =
                            mask_r & (((~(symbol_a)) ^ (symbol_b << rev_counter)) | ((~(left_strand) & top_strand)));
                    rev_combing_cond = ~combing_condition;

                    bitset_left_strand_map[m - 1] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                }

                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);

            }
        }

        //PHASE 1: Process diagonal till fill big left triangle,
        for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

            left_edge = m - 1 - cur_diag_cube_len;
            top_edge = 0;

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

                //parallel  process current
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand >> rev_counter;
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[j];
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


                left_edge++;
                // parallel process previous
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len; ++j) {

                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand << (inside_diag_num + 1);
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[j];
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
                left_edge--;


                //update mask of current move
                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);

            }

//        process last center mask is always all ones
            // parallel process cur
            //todo
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[j])) | ((~left_strand) & top_strand));
                Input rev_combing_cond = ~combing_condition;
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

                //parallel  process current
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < m; ++j) {
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


                top_edge--;
                // parallel process previous
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < m; ++j) {

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
                top_edge++;


                //update mask of current move
                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);

            }




            //proccess center
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
            for (int j = 0; j < m; ++j) {
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) |
                                           ((~left_strand) & top_strand));
                Input rev_combing_cond = ~combing_condition;
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

                //parallel  process current
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len; ++j) {
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


                top_edge--;
                // parallel process previous
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len + 1; ++j) {

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
                top_edge++;


                //update mask of current move
                mask = (mask << 1) | Input(1);
                mask_r = mask_r | (mask_r >> 1);

            }




            //proccess center
#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
            for (int j = 0; j < cur_diag_cube_len; ++j) {
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) |
                                           ((~left_strand) & top_strand));
                Input rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_strand);
                }
            }

        }


        //process last triangle in position  m-1, n-1  cube
        mask = ~Input(0);
        mask_r = mask;
#pragma omp single
        {

            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
                mask = mask << 1;
                mask_r = mask_r >> 1;

                Input left_strand = bitset_left_strand_map[0];
                Input left_cap = left_strand << (inside_diag_num + 1);
                Input symbol_a = a_reverse[0];
                Input symbol_b = b[n - 1];

                Input top_strand = bitset_top_strand_map[n - 1];
                Input combing_condition =
                        mask & ((~((symbol_a << (inside_diag_num + 1)) ^ symbol_b)) | ((~left_cap) & top_strand));
                Input rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_top_strand_map[n - 1] =
                            (rev_combing_cond & top_strand) | (combing_condition & left_cap);

                    top_strand = top_strand >> (inside_diag_num + 1);
                    combing_condition =
                            mask_r &
                            ((~(symbol_a ^ (symbol_b >> (inside_diag_num + 1)))) | ((~(left_strand) & top_strand)));
                    rev_combing_cond = ~combing_condition;
                    bitset_left_strand_map[0] =
                            (rev_combing_cond & left_strand) | (combing_condition & top_strand);

                }


            }
        }

//        aaligned(bitset_top_strand_map,bitset_left_strand_map,a_reverse,b:sizeof(Input)*8) safelen(8

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


//assume both multiple by sizeof(Input)*8
template<class Input>
int prefix_lcs_via_braid_bits_binary(Input *a_reverse, int a_size, int a_total_symbols,
                                     Input *b, int b_size, int b_total_symbols, int alphabet_size) {

    auto bits_per_symbol = int(std::ceil(log2(alphabet_size)));

    Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m) - (m - 1);

    auto not_garbage_bits_a = a_total_symbols % (sizeof(Input) * 8);
    auto not_garbage_bits_b = b_total_symbols % (sizeof(Input) * 8);

    //    at least one bit is okay
    // We set garbage bits of a to  ones and gatbage bits of b to ones  to get d==1 for all garbage bits and then just substract part that always mathes
    if (not_garbage_bits_a) {
        a_reverse[a_size - 1] |= (~((Input(1) << not_garbage_bits_a) - 1));
    }

    if (not_garbage_bits_b) {
        b[b_size - 1] |= (~((Input(1) << not_garbage_bits_b) - 1));
    }


    int left_edge, top_edge;
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

    //PHASE 0:Process first triangle
    for (int inside_diag_num = 0; inside_diag_num <= upper_bound; ++inside_diag_num, rev_counter--) {

        Input left_strand = bitset_left_strand_map[m - 1];
        Input top_strand = bitset_top_strand_map[0];
        Input left_cap = left_strand >> rev_counter;
        Input symbol_a = a_reverse[m - 1];
        Input symbol_b = b[0];
        Input combing_condition =
                mask & (((~(symbol_a >> rev_counter)) ^ symbol_b) | (((~(left_cap)) & top_strand)));
        Input rev_combing_cond = ~combing_condition;

        if (combing_condition) {
            bitset_top_strand_map[0] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            top_strand = top_strand << rev_counter;
            combing_condition =
                    mask_r & (((~(symbol_a)) ^ (symbol_b << rev_counter)) | ((~(left_strand) & top_strand)));
            rev_combing_cond = ~combing_condition;

            bitset_left_strand_map[m - 1] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);
        }

        mask = (mask << 1) | Input(1);
        mask_r = mask_r | (mask_r >> 1);

    }


    //PHASE 1: Process diagonal till fill big left triangle,
    for (int cur_diag_cube_len = 1; cur_diag_cube_len < m; cur_diag_cube_len++) {

        left_edge = m - 1 - cur_diag_cube_len;
        top_edge = 0;

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

            //parallel  process current
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand >> rev_counter;
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[j];
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


            left_edge++;
            // parallel process previous

            for (int j = 0; j < cur_diag_cube_len; ++j) {

                Input left_strand = bitset_left_strand_map[left_edge + j];
                Input top_strand = bitset_top_strand_map[top_edge + j];
                Input left_cap = left_strand << (inside_diag_num + 1);
                Input symbol_a = a_reverse[left_edge + j];
                Input symbol_b = b[j];
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
            left_edge--;


            //update mask of current move
            mask = (mask << 1) | Input(1);
            mask_r = mask_r | (mask_r >> 1);

        }

//        process last center mask is always all ones
        // parallel process cur

        for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[j])) | ((~left_strand) & top_strand));
            Input rev_combing_cond = ~combing_condition;
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
        top_edge = k + 1;//todo

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

            //parallel  process current
            for (int j = 0; j < m; ++j) {
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


            top_edge--;
            // parallel process previous
            for (int j = 0; j < m; ++j) {

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
            top_edge++;


            //update mask of current move
            mask = (mask << 1) | Input(1);
            mask_r = mask_r | (mask_r >> 1);

        }




        //proccess center
        for (int j = 0; j < m; ++j) {
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) |
                                       ((~left_strand) & top_strand));
            Input rev_combing_cond = ~combing_condition;
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

            //parallel  process current
            for (int j = 0; j < cur_diag_cube_len; ++j) {
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


            top_edge--;
            // parallel process previous
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {

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
            top_edge++;


            //update mask of current move
            mask = (mask << 1) | Input(1);
            mask_r = mask_r | (mask_r >> 1);

        }




        //proccess center
        for (int j = 0; j < cur_diag_cube_len; ++j) {
            Input left_strand = bitset_left_strand_map[left_edge + j];
            Input top_strand = bitset_top_strand_map[top_edge + j];
            Input combing_condition = ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) |
                                       ((~left_strand) & top_strand));
            Input rev_combing_cond = ~combing_condition;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) | (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) | (combing_condition & left_strand);
            }
        }

    }


    //process last triangle in position  m-1, n-1  cube
    mask = ~Input(0);
    mask_r = mask;
    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 1;
        mask_r = mask_r >> 1;

        Input left_strand = bitset_left_strand_map[0];
        Input left_cap = left_strand << (inside_diag_num + 1);
        Input symbol_a = a_reverse[0];
        Input symbol_b = b[n - 1];

        Input top_strand = bitset_top_strand_map[n - 1];
        Input combing_condition =
                mask & ((~((symbol_a << (inside_diag_num + 1)) ^ symbol_b)) | ((~left_cap) & top_strand));
        Input rev_combing_cond = ~combing_condition;
        if (combing_condition) {
            bitset_top_strand_map[n - 1] =
                    (rev_combing_cond & top_strand) | (combing_condition & left_cap);

            top_strand = top_strand >> (inside_diag_num + 1);
            combing_condition =
                    mask_r &
                    ((~(symbol_a ^ (symbol_b >> (inside_diag_num + 1)))) | ((~(left_strand) & top_strand)));
            rev_combing_cond = ~combing_condition;
            bitset_left_strand_map[0] =
                    (rev_combing_cond & left_strand) | (combing_condition & top_strand);

        }


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
















//2 bit per symbol

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


//#pragma omp parallel num_threads(threads_num)  default(none) shared(bitset_left_strand_map, bitset_top_strand_map, a_reverse, b, m, n, dis_braid, total_same_length_diag, braid_ones, std::cout)
//    {

        int left_edge, top_edge;
        Input mask;


//#pragma omp  for simd schedule(static) nowait aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < m; ++k) {
            bitset_left_strand_map[k] = braid_ones;
        }

//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map:sizeof(Input)*8)
        for (int k = 0; k < n; ++k) {
            bitset_top_strand_map[k] = Input(0);
        }

        auto upper_bound = (sizeof(Input) * 8 / 2) - 1;
        //process first triangle in 0,0 cube
        mask = Input(1);
        Input mask_r = Input(1) << ((sizeof(Input) * 8) - 2);
        int rev_counter = (sizeof(Input) * 8 - 2);
        //PHASE 0:Process first triangle
//#pragma omp single
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

            left_edge = m - 1 -  cur_diag_cube_len;
            top_edge = 0;
            //to process current
            rev_counter = (sizeof(Input) * 8 - 2);
            mask = Input(1);
            mask_r = Input(1) << rev_counter;

            //to process previous
            Input mask_prev = braid_ones;
            Input mask_prev_r = braid_ones;

//todo

            //process previous size/2 - 1 cubes and current size/2 -1  cubes
            for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num, rev_counter -= 2) {
                //update mask of prev move
                mask_prev <<= 2;
                mask_prev_r >>= 2;

//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
                for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                    Input left_strand = bitset_left_strand_map[left_edge + j];
                    Input top_strand = bitset_top_strand_map[top_edge + j];
                    Input left_cap = left_strand >> rev_counter;
                    Input symbol_a = a_reverse[left_edge + j];
                    Input symbol_b = b[j];


//                    std::cout<<"symbol_a"<<std::bitset<8>(symbol_a)<<std::endl;
//                    std::cout<<"symbol_b"<<std::bitset<8>(symbol_b)<<std::endl;
//                    std::cout<<"symbol_a_>>"<<rev_counter<<"="<<std::bitset<8>(symbol_a>>rev_counter)<<std::endl;


                    Input symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
                    symbols &= (symbols >> 1) & braid_ones;
//                    std::cout<<"symbols"<<std::bitset<8>(symbols)<<std::endl;

                    Input combing_condition =
                            mask & (symbols | (((~(left_cap)) & top_strand)));
                    Input rev_combing_cond = combing_condition ^braid_ones;
//                    std::cout<<"comb"<<std::bitset<8>(combing_condition)<<std::endl;
//                    std::cout<<"Comb reb"<<std::bitset<8>(rev_combing_cond)<<std::endl;



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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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



//    std::cout<<std::bitset<8>(bitset_left_strand_map[0])<<" c "<<std::bitset<8>(bitset_left_strand_map[1])<<std::endl;
//    std::cout<<std::bitset<8>(bitset_top_strand_map[0])<<" c "<<std::bitset<8>(bitset_top_strand_map[1])<<std::endl;
//


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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//#pragma omp  for simd schedule(static) aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//#pragma omp single
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

//        #pragma omp for  simd schedule(static) reduction(+:dis_braid)  aligned(bitset_top_strand_map, bitset_left_strand_map, a_reverse, b:sizeof(Input)*8)
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
//    }

    return a_total_symbols - dis_braid;

}


#endif //CPU_TRANSPOSITION_NETWORK_H
