//
// Created by garrancha on 25.05.2021.
//

#ifndef SEMI_MEMORY_MANAGEMENT_GPU_H
#define SEMI_MEMORY_MANAGEMENT_GPU_H


#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
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

#endif //SEMI_MEMORY_MANAGEMENT_GPU_H
