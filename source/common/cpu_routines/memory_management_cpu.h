//
// Created by garrancha on 25.05.2021.
//

#ifndef SEMI_MEMORY_MANAGEMENT_CPU_H
#define SEMI_MEMORY_MANAGEMENT_CPU_H


template<class Input>
Input *allocate_1D_array_aligned_on_cpu(int num_elements) {
    return static_cast<Input *> (aligned_alloc(sizeof(Input), num_elements * sizeof(Input)));
}

template<class Input>
void free_1D_array_cpu(Input *arr) {
    free(arr);
}

template<class Input>
void copy_from_cpu_to_cpu(Input const *src, Input *dst, int elems) {
    for (int i = 0; i < elems; i++) dst[i] = src[i];
}


#endif //SEMI_MEMORY_MANAGEMENT_CPU_H
