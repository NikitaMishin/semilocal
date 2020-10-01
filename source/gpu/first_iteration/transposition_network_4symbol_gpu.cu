#include <cooperative_groups.h>
#include <cmath>
#include "../memory_management.h"

using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.


template<class Input>
__device__ inline void prefix_braid_init(Input *arr, Input value, int pos) {
    arr[pos] = value;
}


template<class Input>
__device__ inline void process_cube_withoutif(Input &symbol_a, Input &left_strand, Input &symbol_b, Input &top_strand,
                                              Input &l_active_mask, Input &r_active_mask, bool &use_with_mask,
                                              Input &braid_ones) {

    Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;

    int rev_counter = (sizeof(Input) * 8 - 2);
    Input mask = Input(1);

    // upper half
#pragma unroll
    for (int rev_counter = (sizeof(Input) * 8 - 2); rev_counter > 0; rev_counter -= 2) {
        left_cap = left_strand >> rev_counter;
        symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));

        if (use_with_mask) {
            combing_condition &= (l_active_mask >> rev_counter) & r_active_mask;
        }

        rev_combing_cond = combing_condition ^ braid_ones;

        top_strand_shifted = top_strand << rev_counter;
        top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);

        combing_condition <<= rev_counter;
        rev_combing_cond = combing_condition ^ braid_ones;

        left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);

        mask = (mask << 2) | Input(1);
    }


    // center
    symbols = (~(symbol_a ^ symbol_b));
    symbols &= (symbols >> 1) & braid_ones;
    combing_condition = (symbols | ((~left_strand) & top_strand));

    if (use_with_mask) {
        combing_condition &= (l_active_mask) & r_active_mask;
    }


    rev_combing_cond = combing_condition ^ braid_ones;

    top_strand_shifted = top_strand;
    top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_strand);
    left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);

#pragma unroll
    for (int inside_diag_num = 2; inside_diag_num < (sizeof(Input) * 8 / 2 - 1) * 2 + 1; inside_diag_num += 2) {
        mask <<= 2;

        left_cap = left_strand << ((inside_diag_num));
        symbols = ~(((symbol_a << ((inside_diag_num)))) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;

        combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));

        if (use_with_mask) {
            combing_condition &= (l_active_mask << ((inside_diag_num))) & r_active_mask;
        }


        rev_combing_cond = combing_condition ^ braid_ones;

        top_strand_shifted = top_strand >> ((inside_diag_num));
        top_strand = (rev_combing_cond & top_strand) | (combing_condition & left_cap);

        combing_condition >>= ((inside_diag_num));
        rev_combing_cond = combing_condition ^ braid_ones;

        left_strand = (rev_combing_cond & left_strand) | (combing_condition & top_strand_shifted);
    }

}


template<class Input>
__global__ void prefix_braid_withoutif_prestored_lefts(
        Input const *seq_a_rev, int a_size, Input const *seq_b, int b_size,
        Input braid_one,
        Input *bitset_left_strand, Input *bitset_top_strand, int cells_per_thread_l, int cells_per_thread_t,
        int total_thds) {


    int num_diag = a_size + b_size - 1;
    int total_same_length_diag = num_diag - (a_size - 1) - (a_size - 1);
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;

    //sync primitive to sync whole grid
    auto g = this_grid();

    //init_phase
    for (int i = 0; i < cells_per_thread_l; i++) {
        if (total_thds * i + thread_id < a_size)
            prefix_braid_init(bitset_left_strand, braid_one, thread_id + total_thds * i);
    }
    for (int i = 0; i < cells_per_thread_t; i++) {
        if (total_thds * i + thread_id < a_size)
            prefix_braid_init(bitset_top_strand, Input(0), thread_id + total_thds * i);
    }

    Input left_strand = braid_one;
    Input left_symbol = (a_size - 1 - thread_id) >= 0 ? seq_a_rev[thread_id] : 0;
    bool use_with_mask = false;

    g.sync();

    //1 phase
    int b_pos = 0;
    for (int i = a_size - 1; i > 0; i--) {
        // only specific threads  && only active thread should perform
        if (thread_id >= i && thread_id < a_size) {
            Input symbol_b = seq_b[b_pos];
            Input top_strand = bitset_top_strand[b_pos];
            process_cube_withoutif<Input>(left_symbol, left_strand, symbol_b, top_strand, braid_one, braid_one,
                                          use_with_mask,
                                          braid_one);
            b_pos++;
        }
        g.sync();
    }

    //2 phase
    for (int i = 0; i < total_same_length_diag; i++) {
        if (thread_id < a_size) {
            Input symbol_b = seq_b[b_pos];
            Input top_strand = bitset_top_strand[b_pos];
            process_cube_withoutif(left_symbol, left_strand, symbol_b, top_strand, braid_one, braid_one, use_with_mask,
                                   braid_one);
            b_pos++;
        }
        g.sync();
    }

    //3 phase
    for (int i = a_size - 2; i >= 0; i--) {
        if (thread_id <= i) {
            Input symbol_b = seq_b[b_pos];
            Input top_strand = bitset_top_strand[b_pos];
            process_cube_withoutif(left_symbol, left_strand, symbol_b, top_strand, braid_one, braid_one, use_with_mask,
                                   braid_one);
            b_pos++;
        }
        g.sync();
    }

    bitset_left_strand[thread_id] = left_strand;
}


template<class Input>
Input *four_symbol_gpu_runner_fully_gpu(Input *a_reverse_gpu, int a_size, int a_total_symbols,
                                        Input *b_gpu, int b_size, int b_total_symbols,
                                        Input *bitset_left_strand_gpu,
                                        Input *bitset_top_strand_gpu,
                                        int block_size) {
    Input braid_ones = Input(1);

    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
        braid_ones |= (braid_ones << shift * 2);
    }

    int cells_per_thd_l = 1;
    int cells_per_thd_t = std::ceil((1.0 * b_size) / a_size);
    int total_thds = a_size;
    dim3 grid(std::ceil(1.0 * a_size / block_size), 1);
    dim3 block(block_size, 1);

    prefix_braid_withoutif_prestored_lefts <<< grid, block >>>
            (a_reverse_gpu, a_size, b_gpu, b_size, braid_ones, bitset_left_strand_gpu, bitset_top_strand_gpu,
             cells_per_thd_l, cells_per_thd_t, total_thds);
    memory_management::synchronize_with_gpu();

    return bitset_left_strand_gpu;
}







// total thds = a_size




