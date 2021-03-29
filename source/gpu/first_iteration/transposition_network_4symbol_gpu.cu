#include <cooperative_groups.h>
#include <cmath>
#include "../memory_management.h"

using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.


typedef unsigned int Bit32;
typedef unsigned long long Bit64;

namespace bitwise_prefix_lcs_cuda {


    /**
     * Note that we need not to use __syncwarp since we use sync shuffle
     * @param a sum  of  two vectors consisting of 32 32bit values
     * @param b  basically sum of two vectors plus one for each summator
     * @return  Update values for A and B and returns the  upper carry bit  for a value A
     */
    template<bool withIf>
    inline   __device__ int sklansky_sum_reduction(Bit32 &a, Bit32 &b, Bit32 lane_id) {

        Bit32 with_no_carry = a > p; // carry bit of A
        Bit32 with_carry = max(with_no_carry, b > p); // carry bit of B


        Bit32 tmp; // tmp value of A
        Bit32 carry_tmp; // tmp value of carry bit of A
        Bit32 put_other_value; // weather or not we need to place a new value to current variable
        Bit32 active_mask;

        #pragma  unroll
        for (int k = 1; k < 32; k <<= 1) {
            active_mask = (lane_id % k) < k; // threads that on current iteration  can need to swap values

            // tmp values for carries and value
            carry_tmp = with_no_carry;
            tmp = a;

            // update A part
            put_other_value = __shfl_sync(0xffffffff, with_no_carry, k, k << 1); //broadcast
            put_b_value &= active_mask;

            if (withIf) {

                if (put_other_value) {
                    a = b;
                    with_no_carry = with_carry;
                }

            } else {

                a = (a & (put_other_value - 1)) | ((-put_other_value) & b);
                with_no_carry = (~put_other_value_value & with_no_carry) |
                                (put_other_value_value &
                                 with_carry); // if put_b_value == 1 then we put with_no_carry <- with_carry;

            }


            // update B part
            put_other_value = !__shfl_sync(0xffffffff, with_carry, k, k << 1);
            put_other_value &= active_mask;

            if (withIf) {

                if (put_other_value) {
                    b = tmp;
                    with_carry = carry_tmp;
                }

            } else {

                b = (tmp & (put_other_value - 1)) | ((-put_other_value) & b);
                with_carry = (~put_other_value & carry_tmp) | (put_other_value & with_carry);

            }


        }

        return with_no_carry;
    }


    /**
     *
     * @param a_rev
     * @param lookup
     * @param lane_id
     * @param own_carry
     * @param global_id_y
     */
    inline __device__ void kawanami_core(int *a_rev, Bit32 *lookup, int lane_id, Bit32 &own_carry, int global_id_y) {

        int loc = 0; // current position of packed carries for rows
        int bit_pos = 0; // current position of lower bit in loc
        Bit32 carry = 0; // either 1 or 0 for the lower bit of big vector number, for others always 0
        Bit32 carry_pack_to_save = 0; // packed carry bits to save, only upper carry bit is saved
        Bit32 processing_carry_pack = (lane_id == 0) ? own_carry : 0; // only  lower adder can have carry

        #pragma unroll
        for (int i = 0; i < 1024; i++) {

            if (global_id_y < 0 ) {
                // save partial
                // 31st broadcast carry_pack_to_save to all threads
                carry_pack_to_save = __shfl_sync(0xffffffff, carry_pack_to_save, 31, 32);

                if (lane_id == loc) own_carry = carry_pack_to_save; // and loc lane will update own_carry
                return;
            }

            Bit32 key_a = a_rev[global_id_y];
            Bit32 lookup_value = lookup[key_a * n_small + threadIdx.x];
            Bit32 p = vector & lookup_value;

            carry = ((1 << bit_pos) & processing_carry_pack) != 0;  // for others it would be always 0

            Bit32 sum_v = p + vector + carry;
            Bit32 sum_v_inc = sum_v + 1;

            carry_pack_to_save |= (sklansky_sum_reduction<unsigned int, withIf>(sum_v, sum_v_inc) << bit_pos);
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
            global_id_y--;//since a is reverse
        }


    }

    /**
     *
     * possible TODO: add paratrt WIDTH that allows to process more than 32 elements in a row -- cycle optimization
     *
     *
     * @tparam H
     * @tparam W
     * @tparam with_if
     * @param m | Тут все ок, просто offset взять отрицательный
     * @param n_small | 32*Iter
     * @param a
     * @param lookup
     * @param vector_v
     * @param offset_x row offset from left
     * @param offset_y col offset from top
     * @param carries
     */
    template<int H, int BlockDim, int Iter, bool WithIf>
    __global__ void hyrro_kawanami_kernel_without_shared(int *a_rev, int m, Bit32 * lookup, Bit32 *vector_v, Bit32 *carries,
                                                         int n_small, int offset_x, int offset_y) {

        int warp_id = threadIdx.x / 32;
        int lane_id = threadIdx.x % 32;
        int columns_per_block = BlockDim * Iter;// processed by threads within one block since each thread process Iter columns

        // position relative to global
        int global_id_x = offset_x + lane_id + 32 * Iter * warp_id + blockIdx.x * columns_per_block; // global row pos
        int global_id_y = offset_y + H * 1024 * (1 + warp_id + blockIdx.x * (BlockDim / 32)) - 1;// global col pos could be negative to adjust with non divisible parts

        int global_carry_id = global_id_y / 32 - lane_id;

        #pragma unroll
        for (int i = 0; i < Iter; i++) {
            Bit32 vector = vector_v[global_id_x + 32 * i];

            #pragma unroll
            for (int rep = 0; rep < H; rep++) {
                //load packed carries from global memory
                unsigned int own_carry = (global_carry_id - 32 * rep >= 0) ? carries[global_carry_id - 32 * rep] : 0;
                kawanami_core<WithIf>(a_rev, lookup, lane_id, own_carry, global_id_y);

                //save carries for the 1024 elements
                if (global_carry_id - 32 * rep >= 0) carries[global_carry_id] = own_carry;
            }

            // save 1024 bit vector back to memory
            vector_v[global_id_x + 32 * i] = vector;
        }
    }



    template <int Iter, int Split ,int Width, int BlockDim>
    __global__ void hyyro_shared(int *a_rev, int m, Bit32 * lookup, Bit32 *vector_v, Bit32 *carries,
                                 int n_small, int offset_row, int offset_col){

        int per_warp_cols = ((32 - 1) + Width) * Iter;        // per warp processed columns number
        int per_block_cols = per_warp_cols * (BlockDim / 32);     // per block processed columns number


        // Each warp manage own memory space of size WX32. By varying W we make  tradeoffs between occypancy and number of load operations to shared memory
        // to eliminate bank conflict  we fill in column major order
        __shared__ precalced_part Bit32[Split * BlockDim / 32][32];

        //TODO shflup seems to eliminate shared memory here
        __shared__ vector_part Bit32[(BlockDim / 32) *  ((32 - 1) + Split)  ]

        int lane_id = threadIdx.x % 32;
        int warp_id = threadIdx.x / 32;


        //inital stuff
        int global_row = offset_row + lane_id + 32 * Iter * ( warp_id + blockIdx.x * (BlockDim / 32));
        int global_col = offset_col + lane_id + warp_id * per_warp_cols + blockIdx.x * per_block_cols;

        int pos = warp_id * Split;

        #pragma unroll
        for(int k = 0; k < Iter; k++) {

            #pragma unroll
            for(int it = 0; it < Width / Split; it++) {

                // load precalc from shared memory
                #pragma unroll
                for(int i = 0; i < 32; i++) {
                    int key = __shfl_sync(0xffffffff, private_key, i, 32);
                    if (lane_id < Split) precalced_part[lane_id + pos][i] = lookup[key * n_small + global_row + i + it * Split] ;
                }
                // load vectors from global memory
                #pragma unroll
                for (int i = 0; i < 42; i++) {

                    vector_part[lane_id + pos] = vector_v[]



    auto glob_pos = global_id_zero_in_block_col + i * DimBlock;

                    if ((glob_pos >= 0) && (glob_pos < size_b)) {
                        t_part[threadIdx.x + i * DimBlock] = top_strands[glob_pos];
                        b_part[threadIdx.x + i * DimBlock] = b_part[glob_pos];
                    } else {
                        t_part[threadIdx.x + i * DimBlock] = 0;
                        b_part[threadIdx.x + i * DimBlock] = 0;
                    }
                }




                __syncwarp();


                // local vector v

                // jheart
            }



        }


        //load key for each row, no need for shared memory since we can use shfl_sync
        int private_key = (global_row < m) ? a_rev[global_row]: 0;
        Bit32 vector = (global_col < n_small ) ? vector_v[global_col] : 0;



        // load values
        #pragma unroll
        for(int i = 0; i < 32; i++) {
            int key = __shfl_sync(0xffffffff, private_key, i, 32);
            if (lane_id < Split) b_part[ lane_id + pos][i] = lookup[key * n_small + ] ;// if less then 8
        }

    }


}


namespace bitwise_prefix_semi_local_lcs_cuda {



    /**
     * Given 64 unsigned number  0k_1...0k_32 convert it to k_1...k_2  32-bit unsigned int
     * Same idea could be applied to 00k_1....00k_m. Log based approach
     * @param n
     * @return
     */
    __device__ inline unsigned int calc_reduction(unsigned long long n) {

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
     * Cell processing for binary and 4symbol alphabets.
     * l and t both contains 32 strands, so by processing square built on their intersection we process 32*32 cell of initial matrix.
     * Depending on alphabet size, T is either 32 bit unsigned integer (since we can encode each symbol by one bit)
     * or 64 bit unsigned integer (2 bits per symbol).
     * Further, we can apply same technique to higher alphabet sizes, basically uses widen T type and apply more commands to
     * reduce 0^nk_10^nk_t...0 to 32 bit integer k_1...k_32. The changes will be in  calc_reduction.
     * We define parameter K for optimization purposes, the possible values are 0(uint) and 1(ull).
     * It is assumed that there no dummy strands and symbols both in strings and strands. For such specific case see our paper.
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
     * Calculate bitwise semi-local lcs by zig-zag approach.
     * We use additional check of boundaries to allow run kernel even if matrix is small and less then DimGrid
     *
     * Visualization of our zig-zag (stripe) approach. Here we have W = 3 and a warp size of 4. As we can see such approach allows us
     * to each such  quadrangle independently, the only need of sync is needed within quadringle. Latter is achieved by a
     * __syncwarp primitive (on CUDA < 9.0) need not to sync threads within warp at all (what a  good days we left).
     *                  XXX
     *                 XXX
     *                XXX
     *               XXX
     *             ...
     *             ...
     *           XXX
     *          XXX
     *        XXX
     *       XXX
     *    XXX
     *   XXX
     *  XXX
     * XXX
     * @tparam Width how many cells are processed by a single thread
     * @tparam T see  cell_processing definition
     * @tparam K see cell_processing definition
     * @tparam DimBlock same as DimGrid for one Dimensional case (used for ability to use static shared memory allocation)
     * @param a_reversed encoded string, where each symbol requires log2 |Alphabeet| bits, stored in reverse order for consecutive access
     * @param size_a size of input sequence a  in machine words
     * @param b encoded string, where each symbol requires log2 |Alphabeet| bits
     * @param size_b size of input sequence b  in machine words
     * @param left_strands one bit per strand, stored in a reversed order
     * @param top_strands one bit per strand
     * @param offset_b offset from left of the matrix, >= 0
     * @param offset_a offset from bottom of the matrix >=
     */
    template<int Width, class T, int K, int DimBlock>
    __global__ void
    bitwise_prefix_semi_local_kernel(T *a_reversed, int size_a, T *b, int size_b,
                                     unsigned int *left_strands, unsigned int *top_strands, int offset_b,
                                     int offset_a) {

        int per_warp = (32 - 1) + Width;        // per warp processed columns number
        int per_block = per_warp * (DimBlock / 32);    // per block processed columns number


        volatile __shared__ b_part T[per_block]; // we load associated with blocks symbols of b
        volatile __shared__ t_part unsigned int[per_block]; // load part of top strands elements

        auto lane_id = threadIdx.x % 32; // id of thread within his logical warp
        auto warp_id = threadIdx.x / 32; // id of warp that


        auto global_id_col = offset_b + lane_id + warp_id * per_warp + blockIdx.x * per_block;
        auto global_id_zero_in_block_col =
                offset_b + threadIdx.x + blockIdx.x * per_block; // for load to shared memory required stuff

        auto global_id_row = offset_a + threadIdx.x + blockIdx.x * blockDim.x;

        //local registers
        unsigned int l_strand = 0;
        unsigned int t_strand;

        T symbol_a = 0;
        T symbol_b;


        if ((global_id_row >= 0) && (global_id_row < size_a)) {
            symbol_a = a_reversed[global_id_row];
            l_strand = left_strands[global_id_row];
        }

        #pragma unroll
        for (int i = 0; i < (per_block + DimBlock - 1) / DimBlock);
        i++) {
            auto glob_pos = global_id_zero_in_block_col + i * DimBlock;

            if ((glob_pos >= 0) && (glob_pos < size_b)) {
                t_part[threadIdx.x + i * DimBlock] = top_strands[glob_pos];
                b_part[threadIdx.x + i * DimBlock] = b_part[glob_pos];
            } else {
                t_part[threadIdx.x + i * DimBlock] = 0;
                b_part[threadIdx.x + i * DimBlock] = 0;
            }
        }

        int ptr_shared = warp_id * per_warp + lane_id;
        int global_ptr_shared = global_id_col;

        __syncwarp();

#pragma unroll
        for (int i = 0; i < Width; i++) {

            t_strand = t_part[ptr_shared];
            symbol_b = b_part[ptr_shared];

            cell_processing<T, K>(l_strand, t_strand, symbol_a, symbol_b);


            //update t_part
            if ((global_ptr_shared < size_b) && (global_id_row < size_a) && (global_id_row >= 0) &&
                (global_ptr_shared >= 0))
                t_part[ptr_shared] = t_strand;


            __syncwarp();

            ptr_shared++;
            global_ptr_shared++;
        }


        // store l
        if ((global_id_row >= 0) && (lobal_id_row < size_a)) left_strands[global_id_row] = l_strand;


#pragma unroll
        for (int i = 0; i < (per_block + DimBlock - 1) / DimBlock);
        i++) {
            auto glob_pos = global_id_zero_in_block_col + i * DimBlock;
            if ((glob_pos < size_b) && (glob_pos >= 0)) top_strands[glob_pos] = t_part[threadIdx.x + i * DimBlock];
        }

    }




    __device__ inline Bit32 process_hyrro(Bit32 & vector, Bit32 sum_bit , Bit32  lookup_value) {
        Bit32 old_v = vector;
        Bit32 p = lookup_value & old_v;
        with_offset = sum_bit + old_v + p;
        vector =  (old_v ^ p) | with_offset;
        return with_offset < old_v;

    }

}


namespace semi_local_lcs {

    template <int withIf>
    inline  void __device__ process_cell(int symbol_a, int & l, int &t,  int symbol_b) {

        // guess with temporary register will be faster then two xors
        int tmp;

        if (withIf){
            if ((symbol_a == symbol_b) || (l > t) ) {
                tmp = l;
                l = r;
                l = tmp;
            }

        } else {

            bool r = symbol_a == symbol_b || (l > t);
            tmp  = l;

            l = (l & (r - 1))  | ((-r) &  t);
            t = t ^ tmp ^ l;
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




