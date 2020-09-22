#include <cooperative_groups.h>

using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.


/***
 *
 * @param reduced_sticky_braid
 * @param thread_id_global
 * @param offset
 */
__device__ inline void
semi_local_init_phase_conseq_fill_block(int *reduced_sticky_braid, int pos) {
    reduced_sticky_braid[pos] = pos;
}


/**
 *
 * @param reduced_sticky_braid
 * @param seq_a
 * @param seq_b
 * @param left
 * @param top
 * @param left_a
 * @param top_b
 */
__device__ inline void semi_local_fill_phase_diag_without_if(
        int *reduced_sticky_braid, int const *seq_a, int const *seq_b, int left, int top, int left_a, int top_b) {
    // ~xor a b
    int left_strand = reduced_sticky_braid[left_strand];
    int top_strand = reduced_sticky_braid[top_strand];
    int should_swap = (seq_a[left_a] == seq_b[top_b]) || (left_strand > top_strand);
    // aka change
    reduced_sticky_braid[left] = (1 - should_swap) * top_strand + should_swap * left_strand;
    reduced_sticky_braid[top] = (1 - should_swap) * left_strand + should_swap * top_strand;
}


/**
 *
 * @param reduced_sticky_braid
 * @param seq_a
 * @param seq_b
 * @param left
 * @param top
 * @param left_a
 * @param top_b
 */
__device__ inline void semi_local_fill_phase_diag_withif(
        int *reduced_sticky_braid, int const *seq_a, int const *seq_b, int left, int top, int left_a, int top_b) {

    int left_strand = reduced_sticky_braid[left_strand];
    int top_strand = reduced_sticky_braid[top_strand];
    if ((seq_a[left_a] == seq_b[top_b]) || (left_strand > top_strand)) {
        // aka change
        reduced_sticky_braid[left] = top_strand;
        reduced_sticky_braid[top] = left_strand;
    }
}


__global__ void semi_local_gpu_withif(int const *seq_a, int a_size, int const *seq_b, int b_size,
                                      int *reduced_sticky_braid, int cells_per_thread, int total_thds) {
    // a < b
    //thread_id = m


    int num_diag = a_size + b_size - 1;
    int total_same_length_diag = num_diag - (a_size - 1) - (a_size - 1);
    auto thread_id = blockIdx.x * blockDim.x + threadIdx.x;

    //sync primitive to sync whole grid
    auto g = this_grid();

    //init_phase
    for (int i = 0; i < cells_per_thread; i++) {
        if (total_thds * i + thread_id < a_size + b_size)
            semi_local_init_phase_conseq_fill_block(reduced_sticky_braid, thread_id + total_thds * i);
    }
    g.sync();

    for (int i = 0; i < a_size - 1; i++) {

        if (thread_id < a_size + b_size) {
            //todo inline and correct calculation
            if (thread_id >= a_size - 1 - 1 - i) {
                semi_local_fill_phase_diag_withif(reduced_sticky_braid, seq_a, seq_b,
                                                  thread_id, (thread_id + b_size) % (a_size + b_size),
                                                  a_size - 1 - thread_id,
                                                  thread_id % b_size);
            } else {
                semi_local_fill_phase_diag_withif(reduced_sticky_braid, seq_a, seq_b,
                                                  thread_id, thread_id + b_size + i + total_same_length_diag,
                                                  a_size - 1 - thread_id,
                                                  thread_id + i + total_same_length_diag);
            }
        }
        g.sync();
    }

    //phase two
    for (int i = 0; i < total_same_length_diag; i++) {
        semi_local_fill_phase_diag_withif(reduced_sticky_braid, seq_a, seq_b, thread_id,
                                          thread_id + b_size + i, a_size - 1 - thread_id,
                                          thread_id + i);
        g.sync();
    }

}


