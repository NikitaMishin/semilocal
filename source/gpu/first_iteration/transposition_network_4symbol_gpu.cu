#include <cooperative_groups.h>
#include <cmath>
#include "../memory_management.h"

using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.


namespace bitwise_prefix_lcs_cuda {



    /**
     * @tparam T
     * @param a sum  of two sizeof(T)*8 bit unsigned vectors
     * @param b  basically sum + 1
     * @return  Update values for A and B and returns the  upper carry bit  for a value A
     */
    template <class T>
    inline   __device__ int kawanami_sum_reduction_without_if(T &a, T &b) {
        int with_no_carry = a > p;
        int with_carry =  max(with_no_carry, b > p);


        T tmp;
        int carry_tmp;
        int put_other_value;
        int active_mask;

        #pragma  unroll
        for (int k = 1; k < 32; k<<=1) {
            active_mask = (threadIdx.x % k) < k;

            carry_tmp = with_no_carry;
            tmp = a; // save them because  below we assign a possible new values for this variables

            // update A part
            put_other_value = __shfl_up_sync(0xffffffff, with_no_carry, k, k << 1); //broadcast
            put_b_value &= active_mask;


            a = (a & (put_other_value - 1)) | ((-put_other_value) & b);
            with_no_carry = (~put_other_value_value & with_no_carry) |
                            (put_other_value_value & with_carry); // if put_b_value == 1 then we put with_no_carry <- with_carry;

            // update B part
            put_other_value = !__shfl_up_sync(0xffffffff, with_carry, k, k << 1);
            put_other_value &= active_mask ;
            b = (tmp & (put_other_value - 1)) | ((-put_other_value) & b);
            with_carry = (~put_other_value & carry_tmp) | (put_other_value & with_carry);
        }

        return with_no_carry;
    }


    /**
    *
    * @tparam T
    * @param a sum  of two sizeof(T)*8 bit unsigned vectors
    * @param b  basically sum + 1
    * @return  Update values for a and b and returns a upper bit  for a value
    */
    template <class T>
    inline   __device__ int kawanami_sum_reduction_with_if(T &a, T &b) {
        int with_no_carry = a > p;
        int with_carry =  max(with_no_carry, b > p);

        T tmp;
        int carry_tmp;
        int put_other_value;
        int active_mask;

        #pragma  unroll
        for (int k = 1; k < 32; k<<=1) {
            active_mask = (threadIdx.x % k) < k;

            carry_tmp = with_no_carry;
            tmp = a; // save them because  below we assign a possible new values for this variables

            // update A part
            put_other_value = __shfl_up_sync(0xffffffff, with_no_carry, k, k << 1);
            put_other_value &= active_mask ;

            if (put_other_value) {
                a = b;
                with_no_carry = with_carry;
            }

            // update B part
            put_other_value = !__shfl_up_sync(0xffffffff, with_carry, k, k << 1);
            put_other_value &= active_mask;
            if (put_other_value) {
                b = tmp;
                with_carry = carry_tmp;
            }
        }

        return with_no_carry;
    }



    /**
     * With if
     * @tparam T
     * @param m
     * @param n_small | 32
     * @param a
     * @param lookup
     * @param vector_v
     * @param offset_x
     * @param offset_y
     * @param carries | 32
     */
    template <class T>
    __global__ void hyrro_kawanami_kernel(int m, int n_small, int *a  T * lookup, T *vector_v, int offset_x, int offset_y,int *carries) {
// todo volatile
        __shared__ shared_carries T[32];
        __shared__ carry_sum_yes T[32];
        __shared__ carry_sum_no bool[32];
        __shared__ value_A T[32];
        __shared__ value_B T[32];

        //kawanami
        int global_id_x = offset_x + 32 * blockIdx.x + threadIdx.x;
        int global_id_y = 1024 * (offset_y + blockId.x);
        int global_carry_id  = 32 * (offset_y + blockId.x);

        vector = vector_v[global_id_x];


        shared_carries[threadIdx.x] = carries[global_carry_id + threadIdx.x];
        int carry_pack = 0;

        #pragma unroll
        for (int i = 0; i < 1024; i++) {
            //get cur carry

            if(global_id_y > m) break;

            int key_a = a[global_id_y];

            T lookup_value = lookup[key_a * n_small + threadIdx.x];
            T p = vector & lookup_value;
            T sum_v = p + vector +   ;
            T sum_v_inc = sum_v + 1;

            kawanami_sum_reduction()




            value_A[threadIdx.x] = sum_v;
            carry_sum_no[threadIdx.x] = sum_v > p;
            value_B[threadIdx.x] = sum_v + 1;
            carry_sum_yes[threadIdx.x] = max(sum_v > p, (sum_v + 1) > p);
            // reduction

            int n = __shfl_up_sync(0xffffffff, value, i, 8);


            // end todo

            vector = (vector ^ p) | sum_v;
            //carry_pack = todo;

            global_id_y++;
        }

        // save 1024 bit vector back to memory
        vector_v[global_id_x] = vector;
        //save carries for the 1024 elements
        carries[global_carry_id + threadIdx.x] = shared_carries[threadIdx.x];
    }


}


/**
 * Fill array in pos with specified value for each cell
 * @tparam Input
 * @param arr
 * @param value
 * @param pos
 */
template<class Input>
__device__ inline void prefix_braid_init(Input *arr, Input value, int pos) {
    arr[pos] = value;
}


template<class Input>
__device__ inline void process_cube_withoutif(Input &symbol_a, Input &left_strand, Input &symbol_b, Input &top_strand,
                                              Input &l_active_mask, Input &r_active_mask, bool &use_with_mask,
                                              Input &braid_ones) {

    Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;

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
//    auto g = this_grid();
    auto g = this_thread_block();

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
//
    g.sync();

    //1 phase
    int b_pos = 0;
    for (int i = a_size - 1; i > 0; i--) {
        // only specific threads  && only active thread should perform
        if (thread_id >= i && thread_id < a_size) {
            Input symbol_b = seq_b[b_pos];
            Input top_strand = bitset_top_strand[b_pos];
            process_cube_withoutif(left_symbol, left_strand, symbol_b, top_strand, braid_one, braid_one,
                                          use_with_mask,
                                          braid_one);
            bitset_top_strand[b_pos] = top_strand;
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
            bitset_top_strand[b_pos] = top_strand;
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
            bitset_top_strand[b_pos] = top_strand;
            b_pos++;
        }
        g.sync();
    }

    bitset_left_strand[thread_id] = left_strand;
}


template<class Input>
__global__ void
process_diag(
        Input *sticky_braid, Input const *seq_a, int a_size, Input const *seq_b, int b_size, int offset_l,
             Input l_active_mAK,
             int offset_top, int diag_len) {
//
//    int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
//    if (thread_id < diag_len) {
//
//
//        //load left symbol
//        Input a_symbol = seq_a[a_size - offset_l - 1 - thread_id];
//        Input left_strand = seq_a[thread_id + offset_l];
//        Input top_strand = sticky_braid[thread_id + offset_top + a_size];
//        Input b_symbol = seq_b[thread_id + offset_top];
//        bool use_with_mask =
//
//        process_cube_withoutif(a_symbol,left_strand,b_symbol,top_strand,)
//
//
//        sticky_braid[thread_id + offset_l] = (1 - should_swap) * left_strand + should_swap * top_strand;
//        sticky_braid[thread_id + offset_top + a_size] = (1 - should_swap) * top_strand + should_swap * left_strand;
//
//    }


}



template<class Input>
void four_symbol_gpu_runner_fully_gpu(Input *a_reverse_gpu, int a_size, int a_total_symbols,
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
    prefix_braid_withoutif_prestored_lefts <<< grid,block >>>
            (a_reverse_gpu, a_size, b_gpu, b_size, braid_ones, bitset_left_strand_gpu, bitset_top_strand_gpu,
             cells_per_thd_l, cells_per_thd_t, total_thds);
    memory_management::gpuAssert(cudaGetLastError(), __FILE__, __LINE__);

    memory_management::synchronize_with_gpu();

//    return bitset_left_strand_gpu;
}



//
//
//template<class Input>
//Input *four_symbol_fill_gpu_line(Input *a_reverse_gpu, int a_size, int a_total_symbols,
//                                        Input *b_gpu, int b_size, int b_total_symbols,
//                                        Input *bitset_left_strand_gpu,
//                                        Input *bitset_top_strand_gpu,
//                                        int block_size) {
//    Input braid_ones = Input(1);
//
//    for (int shift = 0; shift < sizeof(Input) * 8 / 2; shift++) {
//        braid_ones |= (braid_ones << shift * 2);
//    }
//
//    int cells_per_thd_l = 1;
//    int cells_per_thd_t = std::ceil((1.0 * b_size) / a_size);
//    int total_thds = a_size;
//    dim3 grid(std::ceil(1.0 * a_size / block_size), 1);
//    dim3 block(block_size, 1);
//
//    prefix_braid_withoutif_prestored_lefts<Input> <<< grid, block >>>
//            (a_reverse_gpu, a_size, b_gpu, b_size, braid_ones, bitset_left_strand_gpu, bitset_top_strand_gpu,
//             cells_per_thd_l, cells_per_thd_t, total_thds);
//
//    memory_management::gpuAssert(cudaGetLastError(), __FILE__, __LINE__);
//
//    memory_management::synchronize_with_gpu();
//
//
//    return bitset_left_strand_gpu;
//}
//



// total thds = a_size




