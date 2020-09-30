#include <cooperative_groups.h>

using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.


/***
 * Initial phase. Filling braid:
 * a[i] = i \forall i \in [0...a.size())
 * @param reduced_sticky_braid
 * @param thread_id_global
 * @param offset
 */
__device__ inline void
semi_local_init_phase_conseq_fill(int *reduced_sticky_braid, int pos) {
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
__device__ inline void
semi_local_diag_fill(int *sticky_braid, int &l_strand, int &a_symbol, int const *seq_b, int &top, int &b_pos) {

    int old_left = l_strand;
    int top_strand = sticky_braid[b_pos + top];
    int should_swap = (a_symbol == seq_b[b_pos]) || (l_strand > top_strand);
    // aka change
    l_strand = (1 - should_swap) * l_strand + should_swap * top_strand;
    sticky_braid[b_pos + top] = (1 - should_swap) * top_strand + should_swap * old_left;
}

__global__ void
process_diag(int *sticky_braid, int const *seq_a, int a_size, int const *seq_b, int b_size, int offset_l,
             int offset_top, int diag_len) {

    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
    if (thread_id < diag_len) {


        //load left symbol
        int a_symbol = seq_a[a_size - offset_l - 1 - thread_id];
        int left_strand = seq_a[thread_id + offset_l];
        int top_strand = sticky_braid[thread_id + offset_top + a_size];
        int b_symbol = seq_b[thread_id + offset_top];

        int should_swap = (a_symbol == b_symbol) || (left_strand > top_strand);
        // aka change
        sticky_braid[thread_id + offset_l] = (1 - should_swap) * left_strand + should_swap * top_strand;
        sticky_braid[thread_id + offset_top + a_size] = (1 - should_swap) * top_strand + should_swap * left_strand;

    }


}


__global__ void semi_local_withoutif_prestored_lefts(int const *seq_a, int a_size, int const *seq_b, int b_size,
                                                    int *reduced_sticky_braid, int cells_per_thread, int total_thds) {


    int num_diag = a_size + b_size - 1;
    int total_same_length_diag = num_diag - (a_size - 1) - (a_size - 1);
    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;

    //sync primitive to sync whole grid
    auto g = this_grid();

    //init_phase
    for (int i = 0; i < cells_per_thread; i++) {
        if (total_thds * i + thread_id < a_size + b_size)
            semi_local_init_phase_conseq_fill(reduced_sticky_braid, thread_id + total_thds * i);
    }


    // todo do we need apply forward access pattern and then swap within warp?
    int left_symbol = (a_size - 1 - thread_id) >= 0 ? seq_a[a_size - 1 - thread_id] : 0;
    // economy of
    int left_strand = thread_id;

    g.sync();

    //1 phase
    int b_pos = 0;
    for (int i = a_size - 1; i > 0; i--) {
        // only specific threads  && only active thread should perform
        if (thread_id >= i && thread_id < a_size) {
            semi_local_diag_fill(reduced_sticky_braid, left_strand, left_symbol, seq_b, a_size, b_pos);
            b_pos++;
        }
        g.sync();
    }

    //2 phase
    for (int i = 0; i < total_same_length_diag; i++) {
        if (thread_id < a_size) {
            semi_local_diag_fill(reduced_sticky_braid, left_strand, left_symbol, seq_b, a_size, b_pos);
            b_pos++;
        }
        g.sync();
    }

    //3 phase
    for (int i = a_size - 2; i >= 0; i--) {
        if (thread_id <= i) {
            semi_local_diag_fill(reduced_sticky_braid, left_strand, left_symbol, seq_b, a_size, b_pos);
            b_pos++;
        }
        g.sync();
    }

    reduced_sticky_braid[thread_id] = left_strand;
}




