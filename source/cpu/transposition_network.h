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

//#pragma omp for reduction(+:dis_braid) schedule(static)
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
int prefix_lcs_via_braid_bits_binary(Input *a_reverse, int a_size, int a_total_symbols,
                                     Input *b, int b_size, int b_total_symbols, int alphabet_size) {

    auto bits_per_symbol = int(std::ceil(log2(alphabet_size)));
    auto shift = bits_per_symbol;
    auto word_size_in_bits = sizeof(Input) * 8;

    auto bitset_left_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * a_size));
    auto bitset_top_strand_map = static_cast<Input *> (aligned_alloc(sizeof(Input), sizeof(Input) * b_size));

    auto m = a_size, n = b_size;

    int dis_braid = 0;
    auto num_diag = m + n - 1;
    auto total_same_length_diag = num_diag - (m - 1) - (m - 1);

    auto not_garbage_bits_a = m % (sizeof(Input) * 8);
    auto not_garbage_bits_b = n % (sizeof(Input) * 8);

    //    at least one bit is okay
    // We set garbage bits of a to  ones and gatbage bits of b to ones  to get d==1 for all garbage bits and then just substract part that always mathes
    if (not_garbage_bits_a) {
        a_reverse[a_size - 1] |= (~((Input(1) << not_garbage_bits_a) - 1));
    }

    if (not_garbage_bits_b) {
        b[b_size - 1] |= (~((Input(1) << not_garbage_bits_b) - 1));
    }


    auto threads_num = 1;
//#pragma omp parallel num_threads(threads_num)  default(none) shared(bitset_left_strand_map, bitset_top_strand_map,a_reverse, b,  m, n, dis_braid,total_same_length_diag,std::cout)
//    {

    int left_edge, top_edge;
    Input mask;


//    #pragma omp  for simd schedule(static)
    for (int k = 0; k < m; ++k) {
        bitset_left_strand_map[k] = ~Input(0);
    }

//    #pragma omp  for simd schedule(static)
    for (int k = 0; k < n; ++k) {
        bitset_top_strand_map[k] = Input(0);
    }

    auto upper_bound = (sizeof(Input) * 8) - 1;
    //process first triangle in 0,0 cube
    mask = Input(0);

    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = (mask << 1) | 1;
        auto left_strand = bitset_left_strand_map[m - 1];
        auto top_strand = bitset_top_strand_map[0];
        auto combing_condition = mask & ((~(a_reverse[m - 1] ^ b[0])) | ((~left_strand) & top_strand));
        auto rev_combing_cond = ~combing_condition;
        if (combing_condition) {
            bitset_left_strand_map[m - 1] =
                    (rev_combing_cond & left_strand) || (combing_condition & top_strand);
            bitset_top_strand_map[0] =
                    (rev_combing_cond & top_strand) || (combing_condition & left_strand);
        }
    }

    auto left_strand = bitset_left_strand_map[m - 1];
    auto top_strand = bitset_top_strand_map[0];
    auto combing_condition = ((~(a_reverse[m - 1] ^ b[0])) | ((~left_strand) & top_strand));
    auto rev_combing_cond = ~combing_condition;
    if (combing_condition) {
        bitset_left_strand_map[m - 1] =
                (rev_combing_cond & left_strand) || (combing_condition & top_strand);
        bitset_top_strand_map[0] =
                (rev_combing_cond & top_strand) || (combing_condition & left_strand);
    }


    //    phase one
    for (int cur_diag_cube_len = 0; cur_diag_cube_len < m - 1; cur_diag_cube_len++) {
        mask = Input(0);
        left_edge = m - 1 - cur_diag_cube_len; // todo + 1
        top_edge = 0;
        //        process previous size/2 - 1 cubes and current size/2 -1  cubes with iteration for mask 11111111111
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
            mask = (mask << 1) | 1;

            //parallel  process curent
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                auto left_strand = bitset_left_strand_map[left_edge + j];
                auto top_strand = bitset_top_strand_map[top_edge + j];
                auto combing_condition = mask & ((~(a_reverse[left_edge + j] ^ b[j])) | ((~left_strand) & top_strand));
                auto rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) || (combing_condition & left_strand);
                }
            }
            left_edge--;
            mask = ~mask;
            // parallel process previous
            for (int j = 0; j < cur_diag_cube_len; ++j) {
                auto left_strand = bitset_left_strand_map[left_edge + j];
                auto top_strand = bitset_top_strand_map[top_edge + j];
                auto combing_condition = mask & ((~(a_reverse[left_edge + j] ^ b[j])) | ((~left_strand) & top_strand));
                auto rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) || (combing_condition & left_strand);
                }
            }
            mask = ~mask;
            left_edge++;
        }

//        process last center mask is always all ones
        // parallel process previous
        for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
            auto left_strand = bitset_left_strand_map[left_edge + j];
            auto top_strand = bitset_top_strand_map[top_edge + j];
            auto combing_condition = ((~(a_reverse[left_edge + j] ^ b[j])) | ((~left_strand) & top_strand));
            auto rev_combing_cond = ~combing_condition;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) || (combing_condition & left_strand);
            }
        }
    }


    //phase 2
//    -1 to process  correct
    for (int j = 0; j < total_same_length_diag - 1; ++j) {
        mask = Input(0);
        left_edge = 0;
        top_edge = j;
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
            mask = (mask << 1) | 1;
            mask = ~mask;
            //parallel  process current
            for (int k = 0; k < m; ++k) {
                auto left_strand = bitset_left_strand_map[left_edge + k];
                auto top_strand = bitset_top_strand_map[top_edge + k];
                auto combing_condition =
                        mask & ((~(a_reverse[left_edge + k] ^ b[top_edge + k])) | ((~left_strand) & top_strand));
                auto rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + k] =
                            (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + k] =
                            (rev_combing_cond & top_strand) || (combing_condition & left_strand);
                }
            }
            top_edge++;
            mask =~ mask;
            //parallel  process next
            for (int k = 0; k < m; ++k) {
                auto left_strand = bitset_left_strand_map[left_edge + k];
                auto top_strand = bitset_top_strand_map[top_edge + k];
                auto combing_condition =
                        mask & ((~(a_reverse[left_edge + k] ^ b[top_edge + k])) | ((~left_strand) & top_strand));
                auto rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + k] =
                            (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + k] =
                            (rev_combing_cond & top_strand) || (combing_condition & left_strand);
                }
            }
            top_edge--;

        }
        //last

        //        process last center mask is always all ones
        // parallel process next center
        for (int k = 0; k < m; ++k) {
            auto left_strand = bitset_left_strand_map[left_edge + k];
            auto top_strand = bitset_top_strand_map[top_edge + k];
            auto combing_condition = ((~(a_reverse[left_edge + k] ^ b[top_edge + k])) | ((~left_strand) & top_strand));
            auto rev_combing_cond = ~combing_condition;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + k] =
                        (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + k] =
                        (rev_combing_cond & top_strand) || (combing_condition & left_strand);
            }
        }

    }

    // process cur and prev
    ////    phase 3
    auto start_j = total_same_length_diag;
    for (int cur_diag_cube_len = m - 2; cur_diag_cube_len >= 1; cur_diag_cube_len--, start_j++) {
        mask = Input(0);
        left_edge = 0;
        top_edge = start_j;
        //        process previous size/2 - 1 cubes and current size/2 -1  cubes with iteration for mask 11111111111
        for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
            mask = (mask << 1) | 1;

            //parallel  process curent
            for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
                auto left_strand = bitset_left_strand_map[left_edge + j];
                auto top_strand = bitset_top_strand_map[top_edge + j];
                auto combing_condition =
                        mask & ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) | ((~left_strand) & top_strand));
                auto rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) || (combing_condition & left_strand);
                }
            }
            top_edge--;
            mask = ~mask;
            // parallel process previous
            for (int j = 0; j < cur_diag_cube_len; ++j) {
                auto left_strand = bitset_left_strand_map[left_edge + j];
                auto top_strand = bitset_top_strand_map[top_edge + j];
                auto combing_condition =
                        mask & ((~(a_reverse[left_edge + j] ^ b[top_edge + j])) | ((~left_strand) & top_strand));
                auto rev_combing_cond = ~combing_condition;
                if (combing_condition) {
                    bitset_left_strand_map[left_edge + j] =
                            (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                    bitset_top_strand_map[top_edge + j] =
                            (rev_combing_cond & top_strand) || (combing_condition & left_strand);
                }
            }
            mask = ~mask;
            top_edge++;
        }

//        process last center mask is always all ones
        // parallel process previous
        for (int j = 0; j < cur_diag_cube_len + 1; ++j) {
            auto left_strand = bitset_left_strand_map[left_edge + j];
            auto top_strand = bitset_top_strand_map[top_edge + j];
            auto combing_condition = ((~(a_reverse[left_edge + j] ^ b[j + top_edge])) | ((~left_strand) & top_strand));
            auto rev_combing_cond = ~combing_condition;
            if (combing_condition) {
                bitset_left_strand_map[left_edge + j] =
                        (rev_combing_cond & left_strand) || (combing_condition & top_strand);
                bitset_top_strand_map[top_edge + j] =
                        (rev_combing_cond & top_strand) || (combing_condition & left_strand);
            }
        }
    }



// first process previous level with mask and second process current

    //process first triangle in 0,0 cube
    mask = ~Input(0);
    for (int inside_diag_num = 0; inside_diag_num < upper_bound; ++inside_diag_num) {
        mask = mask << 1;
        auto left_strand = bitset_left_strand_map[0];
        auto top_strand = bitset_top_strand_map[n - 1];
        auto combing_condition = mask & ((~(a_reverse[0] ^ b[n - 1])) | ((~left_strand) & top_strand));
        auto rev_combing_cond = ~combing_condition;
        if (combing_condition) {
            bitset_left_strand_map[0] =
                    (rev_combing_cond & left_strand) || (combing_condition & top_strand);
            bitset_top_strand_map[n - 1] =
                    (rev_combing_cond & top_strand) || (combing_condition & left_strand);
        }
    }


    //TODO is it working only with positive number?
//#pragma omp for simd reduction(+:dis_braid) schedule(static)
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


    //todo think
//    if (std::min(not_garbage_bits_a,not_garbage_bits_b == 0)
    return a_total_symbols - dis_braid;

}
//    return  a_total_symbols  - dis_braid +  (sizeof(Input)*8 - std::max(not_garbage_bits_a,not_garbage_bits_b)) ;


#endif //CPU_TRANSPOSITION_NETWORK_H
