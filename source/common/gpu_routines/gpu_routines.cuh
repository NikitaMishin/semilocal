//
// Created by garrancha on 25.05.2021.
//

#ifndef SEMI_GPU_ROUTINES_CUH
#define SEMI_GPU_ROUTINES_CUH

#include "cpu_routines/memory_management_cpu.h"
#include "gpu_routines/memory_management_gpu.h"

/**
 *
 *  Consider an example: We have  4 threads, a_size < b_size and we process some antidiagonal of size a_size:
 *  _______________________________________
 *  - - - - - - - - ...  - 1 - - - - - - -|
 *  - - - - - - - - ...  0 - - - - - - - -|
 *  - - - - - - - - ...  - - - - - - - - -|
 *  - - - - - - - 2 ...  - - - - - - - - -|
 *  - - - - - - 1 - ...  - - - - - - - - -|
 *  - - - - - 0 - - ...  - - - - - - - - -|
 *  - - - - 3 - - - ...  - - - - - - - - -|
 *  - - - 2 - - - - ...  - - - - - - - - -|
 *  - - 1 - - - - - ...  - - - - - - - - -|
 *  - 0 - - - - - - ...  - - - - - - - - -|
 * @param a
 * @param a_size
 * @param b
 * @param b_size
 * @param raw_kernel
 * @param offset_a
 * @param offset_b
 */

/**
 *
 * @tparam CellsPerThread
 * @param a
 * @param a_size
 * @param b
 * @param b_size
 * @param left_strands
 * @param top_strands
 * @param raw_kernel
 * @param offset_a
 * @param offset_b
 */
template<bool WithIf, bool WithMul>
__global__ void process_antidiagonal_gpu( int *a, int a_size,  int *b, int b_size, int *left_strands,
                                     int *top_strands, int offset_a, int offset_b, int cells_per_thread) {

    int total_threads = blockDim.x * gridDim.x;;
    int global_thread_id = threadIdx.x + blockIdx.x * blockDim.x;
    int row = offset_a + global_thread_id;
    int col = offset_b + global_thread_id;

#pragma unroll
    for (int i = 0; i < cells_per_thread; i++, row += total_threads, col += total_threads) {

        if (row < a_size && col < b_size) {
            int l_strand = left_strands[row];
            int t_strand = top_strands[col];
            int symbol_a = a[row];
            int symbol_b = b[col];
            int cond = (symbol_a == symbol_b) || (l_strand > t_strand);

            if (WithIf) {
                if (cond) {
                    top_strands[col] = t_strand;
                    left_strands[row] = l_strand;
                }
            } else {
                if (WithMul) {
                    // l * cond + t - t * cond = cond * (l  -  t) + t; -- possible overflow for 2^{32-1}
                    top_strands[col] = cond * (l_strand - t_strand) + t_strand;
                    left_strands[row] = cond * (t_strand - l_strand) + l_strand;
                } else {
                    int r_minus = (cond - 1);
                    int minus_r = -cond;
                    left_strands[row] = (l_strand & r_minus) | (minus_r & t_strand);
                    top_strands[col] = (t_strand & r_minus) | (minus_r & l_strand);
                }
            }
        }
    }
}

template<bool WithIf, bool WithMul>
void
process_antidiagonal_wrapper(int  *a, int a_size, int   *b, int b_size, int *left_strands,
                             int *top_strands, int offset_a, int offset_b, int cells_per_thread,
                             int threads_per_block, int diag_len) {

    //setup configuration
    dim3 per_block_thds(std::min(threads_per_block, diag_len), 1, 1);
    dim3 per_grid_blocks(std::ceil(diag_len / per_block_thds.x), 1, 1);
    process_antidiagonal_gpu<WithIf, WithMul><<<per_block_thds, per_grid_blocks>>>(a, a_size, b, b_size,
                                                                               left_strands, top_strands,
                                                                               offset_a, offset_b,cells_per_thread);

}

#endif //SEMI_GPU_ROUTINES_CUH
