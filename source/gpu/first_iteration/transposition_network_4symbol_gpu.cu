#include <cooperative_groups.h>
#include <cmath>
#include "../memory_management.h"

using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.


namespace bitwise_prefix_lcs_cuda {


    /**
     * @tparam T possible types are uint32 and uint64
     * @param a sum  of  two vectors
     * @param b  basically sum of two vectors plus one for each summator
     * @return  Update values for A and B and returns the  upper carry bit  for a value A
     */
    template<class T>
    inline   __device__ int kawanami_sum_reduction_without_if(T &a, T &b) {

        int with_no_carry = a > p; // carry bit of A
        int with_carry = max(with_no_carry, b > p); // carry bit of B


        T tmp; // tmp value of A
        int carry_tmp; // tmp value of carry bit of A
        int put_other_value; // weather or not we need to place a new value to current variable
        int active_mask;

#pragma  unroll
        for (int k = 1; k < 32; k <<= 1) {
            active_mask = (threadIdx.x % k) < k; // threads that on current iteration  can need to swap values

            // tmp values for carries and value
            carry_tmp = with_no_carry;
            tmp = a;

            // update A part
            put_other_value = __shfl_sync(0xffffffff, with_no_carry, k, k << 1); //broadcast
            put_b_value &= active_mask;


            a = (a & (put_other_value - 1)) | ((-put_other_value) & b);
            with_no_carry = (~put_other_value_value & with_no_carry) |
                            (put_other_value_value &
                             with_carry); // if put_b_value == 1 then we put with_no_carry <- with_carry;

            // update B part
            put_other_value = !__shfl_sync(0xffffffff, with_carry, k, k << 1);
            put_other_value &= active_mask;
            b = (tmp & (put_other_value - 1)) | ((-put_other_value) & b);
            with_carry = (~put_other_value & carry_tmp) | (put_other_value & with_carry);
        }

        return with_no_carry;
    }


    template<class T>
    inline   __device__ int kawanami_sum_reduction_with_if(T &a, T &b) {
        int with_no_carry = a > p;
        int with_carry = max(with_no_carry, b > p);

        T tmp;
        int carry_tmp;
        int put_other_value;
        int active_mask;

        #pragma  unroll
        for (int k = 1; k < 32; k <<= 1) {
            active_mask = (threadIdx.x % k) < k;

            carry_tmp = with_no_carry;
            tmp = a;

            // update A part
            put_other_value = __shfl_sync(0xffffffff, with_no_carry, k, k << 1);
            put_other_value &= active_mask;

            if (put_other_value) {
                a = b;
                with_no_carry = with_carry;
            }

            // update B part
            put_other_value = !__shfl_sync(0xffffffff, with_carry, k, k << 1);
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
    __global__ void hyrro_kawanami_kernel_without_shared(int m, int n_small, int *a  unsigned int *
    lookup, unsigned int *vector_v, int offset_x, int offset_y,
                                                         unsigned int *carries) {
        //TODO does it faster then with shared?

        // position relative to global
        int global_id_x = offset_x + 32 * blockIdx.x + threadIdx.x;
        int global_id_y = 1024 * (offset_y + blockId.x);
        int global_carry_id = 32 * (offset_y + blockId.x);


        //load packed carries from global memory
        unsigned int own_carry = carries[global_carry_id + threadIdx.x];

        unsigned int vector = vector_v[global_id_x];


        int loc = 0; // current position of packed carries for rows
        int bit_pos = 0; // current position of lower bit in loc
        unsigned int carry = 0; // either 1 or 0 for the lower bit of big vector number, for others always 0
        unsigned int carry_pack_to_save = 0; // packed carry bits to save, only upper carry bit is saved
        unsigned int processing_carry_pack = (threadIdx.x == 0) ? shared_carries[loc]
                                                                : 0; // only  lower adder can have carry

        #pragma unroll
        for (int i = 0; i < 1024; i++) {

            if (global_id_y > m) {
                // save partial
                // 31st broadcast carry_pack_to_save to all threads
                carry_pack_to_save = __shfl_sync(0xffffffff, carry_pack_to_save, 31, 32);
                if (threadIdx.x == loc) own_carry = carry_pack_to_save; // and loc lane will update own_carry

                break; // out of matrix
            }

            int key_a = a[global_id_y];
            T lookup_value = lookup[key_a * n_small + threadIdx.x];
            T p = vector & lookup_value;


            carry = ((1 << bit_pos) & processing_carry_pack) != 0;  // for others it would be always 0

            T sum_v = p + vector + carry;
            T sum_v_inc = sum_v + 1;

            carry_pack_to_save |= (kawanami_sum_reduction_without_if(sum_v, sum_v_inc) << bit_pos);
            vector = (vector ^ p) | sum_v; // update vector

            // if 31 then we need to save packed values and load another one
            if ((i % 32) == 31) {

                // 31st broadcast carry_pack_to_save to all threads
                carry_pack_to_save = __shfl_sync(0xffffffff, carry_pack_to_save, 31, 32);
                if (threadIdx.x == loc) own_carry = carry_pack_to_save; // and loc lane will update own_carry
                carry_pack_to_save = 0;

                loc++;
                // transfer carry pack from loc thread to  lane with id = loc
                processing_carry_pack = __shfl_sync(0xffffffff, own_carry, loc, 32);
                processing_carry_pack = (threadIdx.x == 0 && loc < 32) ? processing_carry_pack : 0;

                bit_pos = -1; // will be 0 at the end of this loop iteration
            }

            bit_pos++;
            global_id_y++;
        }

        // save 1024 bit vector back to memory
        vector_v[global_id_x] = vector;
        //save carries for the 1024 elements
        carries[global_carry_id + threadIdx.x] = own_carry;

    }


}


namespace bitwise_prefix_semi_local_lcs_cuda {


    /**
     * Given 64 unsigned number  0k_1...0k_32 convert it to k_1...k_2  32-bit unsigned int
     * Same idea could be applied to 00k_1....00k_m. Log based approach
     * @param n
     * @return
     */
    __device__ inline unsigned int calc_reduction(unsigned long long int n) {

        n &= (n >> 1) & (6148914691236517205ull);

        n |= (n >> 1);
        n &= (3689348814741910323ull);
        n |= ((n >> 2));
        n &= (1085102592571150095ull);
        n |= ((n >> 4));
        n &= (71777214294589695ull);
        n |= ((n >> 8));
        n &= (281470681808895ull);
        n |= ((n >> 16));
        n &= (4294967295ull);

        return unsigned int(n);
    }


    /**
     * Id function
     * @param n
     * @return
     */
    __device__ inline unsigned int calc_reduction(unsigned int n) {
        return n;
    }


    /**
     * Cell processing for binary and 4symbol
     *
     * @tparam T possible types are unsigned int (for binary), ull  for 4 symbol alphabet; for 16 symbol see paper
     * @tparam K possible values are 0 (unsigned int), 1 (unsigned long long)
     * @param l left packed 32 strands in machine word
     * @param t top packed 32 strands in machine word
     * @param a packed symbols of string a of type T (holds 32 symbols)
     * @param b packed symbols of string b of type T (holds 32 symbols)
     * @return updated values of l and t
     */
    template<class T, int K>
    __device__ inline void cell_processing(unsigned int &l, unsigned int &t, T a, T b) {

        unsigned int l_shifted = l;
        unsigned int t_shifted = t;
        unsigned int mask = 1;

        unsigned int cond;


#pragma unroll
        for (int shift = 31; shift > 0; shift--) {

            l_shifted = l >> shift;
            t_shifted = t << shift;

            cond = calc_reduction(~((a >> (shift << K)) ^ b));

            t = (l_shifted | (~mask)) & (t | (cond & mask));
            l = t_shifted ^ (t << shift) ^ l;

            mask = (mask << 1) | Input(1);
        }

        cond = calc_reduction(~(a ^ b));

        l_shifted = l;
        t_shifted = t;

        t = (l_shifted | (~mask)) & (t | (cond & mask));
        l = t_shifted ^ (t) ^ l;

        mask = unsigned int(-1);

#pragma unroll
        for (int shift = 1; shift < 32; shift++) {
            mask <<= 1;

            l_shifted = l << shift;
            t_shifted = t >> shift;

            cond = calc_reduction(~(((a << (shift << K)) ^ b)));

            t = (l_shifted | (~mask)) & (t | (cond & mask));
            l = t_shifted ^ (t >> shift) ^ l;
        }
    }



    /**
     *
     * @tparam Width  number of processing cells per thread. Witdh of one is a special case of antidiagonal patttern
     * @tparam T
     * @tparam K
     * @param a_reversed
     * @param size_a
     * @param b
     * @param size_b
     * @param left_strands
     * @param top_strands
     * @param offset_b
     * @param offset_a
     */
    template<int Width, class T, int K, DimBlock>
    __global__ void
    bitwise_prefix_semi_local_kernel(T *a_reversed, int size_a, T *b, int size_b,
                  unsigned int *left_strands, unsigned int *top_strands, int offset_b, int offset_a) {

        int per_warp = (32 - 1) + Width;        // per warp
        int per_block = per_warp * DimBlock;    // per block processed


        //todo dynamic shared
        volatile __shared__ b_part T[per_block];
        volatile __shared__ t_part unsigned int[per_block];

        auto lane_id = threadIdx.x % 32;
        auto warp_id = threadIdx.x / 32;


        auto global_id_col = offset_a +  lane_id + warp_id  * per_warp +   blockIdx.x * per_block;

        auto global_id_zero_in_block_col = offset_a + blockIdx.x * per_block + threadId.x;
        auto global_id_row = offset_b + threadIdx.x + blockIdx.x * blockDim.x;

        unsigned int l_strand = 0;
        unsigned int t_strand;

        T symbol_a = 0;
        T symbol_b;


        if (global_id_row >= 0 & global_id_row < size_a) {
            symbol_a = a_reversed[global_id_row];
            l_strand = left_strands[global_id_row];
        }

        #pragma unroll
        for (int i = 0; i < (per_block + DimBlock - 1) / DimBlock); i++) {
            auto glob_pos = global_id_zero_in_block_col + i * DimBlock;
            // todo i guess it is always >= 0 or not
            if ( (glob_pos >= 0) & ( glob_pos < size_b) ) {
                t_part[threadIdx.x+ i * DimBlock] = top_strands[glob_pos];
                b_part[threadIdx.x+ i * DimBlock] = b_part[glob_pos];
            }
        }

        int ptr_shared = warp_id * per_warp + lane_id;
        int global_ptr_shared = global_id_col;


        #pragma unroll
        for (int i = 0; i < Width; i++) {

            // todo is last condition really necessary?
            if (global_ptr_shared > 0 & global_ptr_shared < size_b & global_id_row > 0) {

                t_strand = t_part[ptr_shared];
                symbol_b = b_part[ptr_shared];

                cell_processing<T, K>(l_strand, t_strand, symbol_a, symbol_b);

                //update t_part
                t_part[ptr_shared] = t_part;
            }

            ptr_shared++;
            global_ptr_shared++;
        }


        // store l
        if (global_id_row >= 0 & global_id_row < size_a) left_strands[global_id_row] = l_strand;


        #pragma unroll
        for (int i = 0; i < (per_block + DimBlock - 1) / DimBlock); i++) {
            auto glob_pos = global_id_zero_in_block_col + i * DimBlock;
            if ( (glob_pos >= 0) & ( glob_pos < size_b) ) top_strands[glob_pos] t_part[threadIdx.x+ i * DimBlock];
        }

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
    prefix_braid_withoutif_prestored_lefts <<< grid, block >>>
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




