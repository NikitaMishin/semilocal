//
// Created by garrancha on 17.05.2021.
//

#ifndef GPU_ALGORITHMS_CUH
#define GPU_ALGORITHMS_CUH

#include "utils/types.h"

#include "monge/matrices.h"


namespace semi_local {
    namespace gpu {


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
        __global__ void process_antidiagonal(const int *a, int a_size, const int *b, int b_size, int *left_strands,
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
        process_antidiagonal_wrapper(const GpuInt *a, int a_size, const GpuInt *b, int b_size, GpuInt *left_strands,
                                     GpuInt *top_strands, int offset_a, int offset_b, int cells_per_thread,
                                     int threads_per_block, int diag_len) {

            //setup configuration
            dim3 per_block_thds(std::min(threads_per_block, diag_len), 1, 1);
            dim3 per_grid_blocks(std::ceil(diag_len / per_block_thds.x), 1, 1);
            process_antidiagonal<WithIf, WithMul><<<per_block_thds, per_grid_blocks>>>(a, a_size, b, b_size,
                                                                                       left_strands, top_strands,
                                                                                       offset_a, offset_b);

        }

    }


    namespace implementation {


    }


}


template<class Input>
class SemiLocalStrategy {
public:
    virtual void compute(AbstractPermutation &permutation, const Input *a, int a_size, const Input *b, int b_size) = 0;
};

template<class Input>
class GPUSemiLocalStrategy : public SemiLocalStrategy<Input> {
};

template<class Input, bool WithIf, bool WithMul>
class GpuSemiLocalAntiDiagonal : public GPUSemiLocalStrategy<Input> {
private:
    int thread_per_block;

public:
    explicit GpuSemiLocalAntiDiagonal(int thds_per_block) : thread_per_block(thds_per_block) {}


    void compute(AbstractPermutation &permutation, const Input *a, int a_size, const Input *b, int b_size) override {
        using namespace semi_local::gpu;
        using namespace memory_management;

        auto host_a = allocate_1D_array_aligned_on_cpu<Input>(a_size);
        copy_from_cpu_to_cpu(a, host_a, a_size);
        auto host_b = allocate_1D_array_aligned_on_cpu<Input>(b_size);
        copy_from_cpu_to_cpu(b, host_b, b_size);

        auto host_l_strands = allocate_1D_array_aligned_on_cpu<Input>(a_size);
        auto host_t_strands = allocate_1D_array_aligned_on_cpu<Input>(b_size);

        for (int i = 0; i < a_size; i++) host_l_strands[i] = i;
        for (int i = 0; i < b_size; i++) host_t_strands[i] = a_size + i;


        auto device_a = allocate_1D_array_aligned_on_gpu<Input>(a_size);
        auto device_b = allocate_1D_array_aligned_on_gpu<Input>(b_size);
        auto device_l_strands = allocate_1D_array_aligned_on_cpu<Input>(a_size);
        auto device_t_strands = allocate_1D_array_aligned_on_cpu<Input>(b_size);

        copy_from_cpu_to_gpu_sync(host_a, device_a, a_size);
        copy_from_cpu_to_gpu_sync(host_b, device_b, b_size);
        copy_from_cpu_to_gpu_sync(host_l_strands, device_l_strands, a_size);
        copy_from_cpu_to_gpu_sync(host_t_strands, device_t_strands, b_size);


        // phase 1
        auto size_len = 1;
        int offset_a, offset_b, cells_per_thread;
        for (int i = 0; i < std::min(a_size, b_size) - 1; i++, size_len++) {
            process_antidiagonal_wrapper<WithIf, WithMul>(device_a, a_size, device_b, b_size, device_l_strands,
                                                          device_t_strands, offset_a, offset_b, cells_per_thread,
                                                          thread_per_block, size_len);
        }


    }

};


#endif //GPU_ALGORITHMS_CUH
