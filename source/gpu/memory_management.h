//
// Created by nikita on 01.10.2020.
//

#ifndef MEMORY_MANAGEMENT_H
#define MEMORY_MANAGEMENT_H


#include <cstdio>
#include <iostream>


namespace memory_management {

    #define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }
    inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
        if (code != cudaSuccess) {
            fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
            if (abort) exit(code);
        }
    }

    template<class Input>
    Input *allocate_1D_array_aligned_on_cpu(int num_elements) {
        return static_cast<Input *> (aligned_alloc(sizeof(Input), num_elements * sizeof(Input)));
    }

    template<class Input>
    void  free_1D_array_cpu(Input *arr) {
        free(arr);
    }


    template<class Input>
    Input *allocate_1D_array_aligned_on_gpu(int num_elements) {
        Input *arr;
        gpuErrchk(cudaMalloc((void **) &arr, sizeof(Input) * num_elements));
        return arr;
    }

    template<class Input>
    void free_1D_array_gpu(Input *arr) {
        gpuErrchk(cudaFree(arr));
    }

    template<class Input>
    void copy_from_cpu_to_gpu_sync(Input *src_cpu, Input *dst_gpu, int elems) {
        gpuErrchk(cudaMemcpy(dst_gpu, src_cpu, sizeof(Input) * elems, cudaMemcpyHostToDevice));
    }

    template<class Input>
    void copy_from_gpu_to_cpu_sync(Input *src_gpu, Input *dst_cpu, int elems) {
        gpuErrchk(cudaMemcpy(dst_cpu, src_gpu, sizeof(Input) * elems, cudaMemcpyDeviceToHost));
    }

    void synchronize_with_gpu() {
        cudaDeviceSynchronize();
    }
    template<class Input>
    void copy_from_cpu_to_cpu(Input * src, Input * dst, int elems) {
        for (int i = 0; i < elems; i++) dst[i] = src[i];
    }


}


#endif //MEMORY_MANAGEMENT_H
