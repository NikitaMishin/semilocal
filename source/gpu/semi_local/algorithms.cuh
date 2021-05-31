//
// Created by garrancha on 17.05.2021.
//

#ifndef GPU_ALGORITHMS_CUH
#define GPU_ALGORITHMS_CUH

#include "utils/types.h"

#include "monge/matrices.h"
#include "cpu_routines/memory_management_cpu.h"
#include "gpu_routines/memory_management_gpu.h"
#include "cpu_routines/cpu_routines.h"
#include "gpu_routines/gpu_routines.cuh"


namespace semi_local {


}


template<class Input>
class SemiLocalStrategy {
public:
    virtual void compute(AbstractPermutation &permutation,  Input const *a, int a_size,  Input const *b, int b_size) = 0;
};

template<class Input>
class GPUSemiLocalStrategy : public SemiLocalStrategy<Input> {
};



/**
 *
 * @tparam Input
 * @tparam WithIf
 * @tparam WithMul
 */
template<class Input,bool GPUEnable, bool WithIf, bool WithMul>
class SemiLocalAntiDiagonalStrategy : public GPUSemiLocalStrategy<Input> {
private:
    int thds_per_block_gpu;
    int cpu_lower_bound;
    int cpu_threads_count;
    int cells_per_thds_threshold;

public:

    explicit SemiLocalAntiDiagonalStrategy(int cpu_threads_count,int threads_per_block_gpu = 512,
                                           int len_to_switch_to_gpu = 1000, int cells_per_thds_gpu_threshold = 15000): cpu_threads_count(cpu_threads_count),
                                                                             thds_per_block_gpu(threads_per_block_gpu),
                                                                             cpu_lower_bound(len_to_switch_to_gpu),
                                                                             cells_per_thds_threshold(cells_per_thds_gpu_threshold){}



    void compute(AbstractPermutation &permutation,  Input const * a, int a_size, Input  const*b, int b_size) override {
        if (a_size > b_size){
            auto tmp = Permutation(a_size+b_size,a_size+b_size);
            compute(tmp,b,b_size,a,a_size);
            fill_permutation_ba(&tmp, &permutation, a_size, b_size);
            return;
        }
        // now a_size <= b_size

        auto host_l_strands = allocate_1D_array_aligned_on_cpu<Input>(a_size);
        auto host_t_strands = allocate_1D_array_aligned_on_cpu<Input>(b_size);

        auto host_b = allocate_1D_array_aligned_on_cpu<Input>(b_size);
        copy_from_cpu_to_cpu(b, host_b, b_size);
        auto host_a_rev = allocate_1D_array_aligned_on_cpu<Input>(a_size);
        #pragma omp parallel num_threads(cpu_threads_count)
        {
            fill_a_reverse_cpu(a,host_a_rev,a_size);
            initialization_cpu(host_l_strands, host_t_strands, a_size, b_size);
        }


        #pragma omp parallel num_threads(cpu_threads_count)
        for (int i = 1; i < std::min(a_size, cpu_lower_bound); i++) {
            anti_diagonal_computation_cpu<Input,WithIf,true>(host_l_strands, host_t_strands, host_a_rev, host_b,
                                                                 i, a_size - i, 0);
        }

        Input * device_a_rev;
        Input * device_b;
        Input * device_l_strands;
        Input * device_t_strands;

        // we switch to gpu computation if GPU enabled and lower bound is hit
        bool should_switch = a_size >= cpu_lower_bound;

        if (GPUEnable && should_switch) {
            device_a_rev = allocate_1D_array_aligned_on_gpu<Input>(a_size);
            device_b = allocate_1D_array_aligned_on_gpu<Input>(b_size);
            device_l_strands = allocate_1D_array_aligned_on_gpu<Input>(a_size);
            device_t_strands = allocate_1D_array_aligned_on_gpu<Input>(b_size);

            copy_from_cpu_to_gpu_sync(host_a_rev, device_a_rev, a_size);
            copy_from_cpu_to_gpu_sync(host_b, device_b, b_size);
            copy_from_cpu_to_gpu_sync(host_l_strands, device_l_strands, a_size);
            copy_from_cpu_to_gpu_sync(host_t_strands, device_t_strands, b_size);
        }


        //rest of first phase
        if(!GPUEnable) {
            #pragma omp parallel num_threads(cpu_threads_count)
            for (int i = std::min(a_size, cpu_lower_bound); i < a_size; i++) {
                anti_diagonal_computation_cpu<Input, WithIf, true>(host_l_strands, host_t_strands, host_a_rev, host_b,
                                                                   i, a_size - i, 0);
            }
        } else {
            for (int i = std::min(a_size, cpu_lower_bound); i < a_size; i++) {
                int cell_per_thds = int(ceil(double(i) / cells_per_thds_threshold));
                process_antidiagonal_wrapper<WithIf, WithMul>(device_a_rev, a_size, device_b, b_size, device_l_strands,
                                                              device_t_strands, a_size - i, 0,
                                                              cell_per_thds, thds_per_block_gpu, i);
                synchronize_with_gpu();
            }
        }

        // second phase
        auto second_phase_len = b_size  - (a_size - 1);
        auto len = a_size;
        int cell_per_thds = int(ceil(double(len) / cells_per_thds_threshold));

        if (!GPUEnable) {
            #pragma omp parallel num_threads(cpu_threads_count)
            for (int i = 0; i < second_phase_len; i++) {
                anti_diagonal_computation_cpu<Input, WithIf, true>(host_l_strands, host_t_strands, host_a_rev, host_b,
                                                                   len, 0, i);
            }
        } else {
            for (int i = 0; i < second_phase_len; i++) {
                process_antidiagonal_wrapper<WithIf, WithMul>(device_a_rev, a_size, device_b, b_size, device_l_strands,
                                                              device_t_strands, 0, i,
                                                              cell_per_thds, thds_per_block_gpu, len);
                synchronize_with_gpu();
            }
        }

        int offset = second_phase_len;
        //third phase
        if (!GPUEnable) {
            #pragma omp parallel num_threads(cpu_threads_count)
            for (int diag_len = a_size - 2; diag_len >= 0; --diag_len) {
                anti_diagonal_computation_cpu<Input, WithIf, true>(host_l_strands,host_t_strands,host_a_rev,host_b,
                                                                   diag_len+1,0, offset);
                offset++;
            }
        } else {
            for (int diag_len = a_size - 2; diag_len > cpu_lower_bound; --diag_len ) {
                process_antidiagonal_wrapper<WithIf, WithMul>(device_a_rev, a_size, device_b, b_size, device_l_strands,
                                                              device_t_strands, 0, offset, cell_per_thds,
                                                              thds_per_block_gpu, diag_len + 1);
                synchronize_with_gpu();
                offset++;
            }

            if (GPUEnable && should_switch) {
                copy_from_gpu_to_cpu_sync(device_t_strands,host_t_strands,b_size);
                copy_from_gpu_to_cpu_sync(device_l_strands,host_l_strands,a_size);

                free_1D_array_gpu(device_l_strands);
                free_1D_array_gpu(device_t_strands);
                free_1D_array_gpu(device_b);
                free_1D_array_gpu(device_a_rev);
            }

            #pragma omp parallel num_threads(cpu_threads_count)
            for (int diag_len = cpu_lower_bound; diag_len >= 0; --diag_len ) {
                anti_diagonal_computation_cpu<Input, WithIf, true>(host_l_strands,host_t_strands,host_a_rev,host_b,
                                                                   diag_len+1,0, offset);
                offset++;
            }
        }
        //end of third phase


        #pragma omp parallel num_threads(cpu_threads_count)
        {
            construct_permutation_cpu(permutation, host_l_strands, host_t_strands, false, a_size, b_size);
        }

        free_1D_array_cpu(host_a_rev);
        free_1D_array_cpu(host_b);
        free_1D_array_cpu(host_l_strands);
        free_1D_array_cpu(host_t_strands);

    }
};


#endif //GPU_ALGORITHMS_CUH
