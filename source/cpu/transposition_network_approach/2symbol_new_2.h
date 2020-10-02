//
// Created by nikita on 02.10.2020.
//

#ifndef CPU_2SYMBOL_NEW_2_H
#define CPU_2SYMBOL_NEW_2_H


#include <vector>
#include <cmath>
#include <bitset>

template<class Input, class SymbolType>
inline Input calc_reduction(SymbolType & n){
    n |= (n >> 1);
    n &= (3689348814741910323ll);
    n |= ((n >> 2));
    n &= (1085102592571150095ll);
    n |= ((n >> 4));
    n &= (71777214294589695ull);
    n |= ((n >> 8));
    n &= (281470681808895ll);
    n |= ((n >> 16));
    n &= (4294967295ll);
    return Input(n);
}

template<class Input, class SymbolType>
inline void process_cubes_antidiag_mpi_bin_types_two(int lower_bound, int upper_bound, int left_edge, int top_edge,
                                           Input *bitset_left_strand_map,
                                           Input *bitset_top_strand_map,
                                           SymbolType *a_reverse, SymbolType *b) {

    const int upper = sizeof(Input) * 8 - 1;

#pragma omp   for  simd schedule(static)  aligned(bitset_top_strand_map, bitset_left_strand_map:sizeof(Input)*8) aligned(a_reverse, b:sizeof(SymbolType)*8)
    for (int j = lower_bound; j < upper_bound; ++j) {

        Input left_cap, combing_condition, rev_combing_cond, top_strand_shifted;
        Input left_strand = bitset_left_strand_map[left_edge + j];
        Input top_strand = bitset_top_strand_map[top_edge + j];

        SymbolType symbol_a = a_reverse[left_edge + j];
        SymbolType symbol_b = b[top_edge + j];

        SymbolType symbols;
        Input res_symbols;

        Input mask = Input(1);


        // upper half
#pragma GCC unroll  256
        for (int rev_counter = (sizeof(Input) * 8 - 1); rev_counter > 0; rev_counter--) {
            left_cap = left_strand >> rev_counter;
            symbols = ~(((symbol_a >> (rev_counter * 2))) ^ symbol_b);

            symbols &= (symbols >> 1) & (6148914691236517205ull);
            res_symbols = calc_reduction<Input, SymbolType>(symbols);

            combing_condition = mask & (res_symbols | (((~(left_cap)) & top_strand)));
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
        symbols &= (symbols >> 1) & (6148914691236517205ULL);
        res_symbols = calc_reduction<Input, SymbolType>(symbols);


        combing_condition = (res_symbols | ((~left_strand) & top_strand));
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
            symbols = ~(((symbol_a << ((inside_diag_num * 2))) ^ symbol_b));
            symbols &= (symbols >> 1) & (6148914691236517205ULL);

            res_symbols = calc_reduction<Input, SymbolType>(symbols);

            combing_condition = mask & (res_symbols | (((~(left_cap)) & top_strand)));
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


template<class Input, class SymbolType>
int prefix_lcs_via_braid_bits_2symbol_v3_full_mask(SymbolType *a_reverse, int a_size, int a_total_symbols,
                                                   SymbolType *b, int b_size, int b_total_symbols, int threads_num) {


     Input *bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    Input *bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));


    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);



#pragma omp parallel num_threads(threads_num) default(none) shared(bitset_left_strand_map, bitset_top_strand_map, a_reverse, b, m, n, dis_braid, total_same_length_diag)
    {

#pragma omp  for simd schedule(static)  aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < n; ++k) {
            bitset_top_strand_map[k] = Input(0);
        }

#pragma omp  for simd schedule(static) aligned(bitset_left_strand_map:sizeof(Input)*8)
        for (int k = 0; k < m; ++k) {
            bitset_left_strand_map[k] = ~Input(0);
        }

        for (int diag_len = 0; diag_len < m - 1; diag_len++) {
            process_cubes_antidiag_mpi_bin_types_two<Input,SymbolType>(0, diag_len + 1, m - 1 - diag_len, 0, bitset_left_strand_map,
                                           bitset_top_strand_map, a_reverse, b);

        }

        for (int k = 0; k < total_same_length_diag; k++) {
            process_cubes_antidiag_mpi_bin_types_two<Input,SymbolType>(0, m, 0, k, bitset_left_strand_map,
                                           bitset_top_strand_map, a_reverse, b);
        }

        auto start_j = total_same_length_diag;

        for (int diag_len = m - 1; diag_len >= 1; diag_len--) {
            process_cubes_antidiag_mpi_bin_types_two<Input,SymbolType>(0, diag_len, 0, start_j, bitset_left_strand_map,
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


#endif //CPU_2SYMBOL_NEW_2_H
