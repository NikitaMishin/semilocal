//
// Created by garrancha on 11.05.2021.
//

#ifndef CPU_KAWANAMI_H
#define CPU_KAWANAMI_H

#include "utils/types.h"

namespace bitwise_prefix_lcs_cuda {

//
// <---->
//
//
//



    namespace hyyro {


        namespace stripe {

            /**
             *
             * Example ChunkSize = 3, Width = 6, warp size = 3, Iter = 2
             *
             *                             ___
             *                             XXXXXX|        | a_size ....                             b_size
             *                            XXXXXX |        |
             *                           XXXXXX  |        |
             *                     XXXXXX        |        |
             *                    XXXXXX         |<- warp |
             *                   XXXXXX          |        |
             *                                            |
             *                                            |
             *                                            |
             *                                            |
             *                                            |
             *                                            |
             *                                            |  <- threadblock
             *          |columns |                        |
             *          | per    |                        |
             *          | warp   |                        |
             *          |------- |                        |
             *          |  width |                        |
             *  ...     |  ----- |                        |
             *          |  XXXXXX|            |           |
             *          | XXXXXX |            |           |
             *          |XXXXXX  |            |Iter * warp|
             *     XXXXXX        |            |           |
             *    XXXXXX         |<-warp_0    |           |
             *   XXXXXX          |            |           | 0
             *
             *
             *
             *
             * @tparam Iter number of Iterations for each warp
             * @tparam ChunkSize  Width is divisile by ChunkSize. We split one iteration on Width/ChunkSize for reduce shared memory reqieremetns to allow better occupancy
             * @tparam Width number of columns that  is processed by each thread in one Iter
             * @tparam BlockDim same as BlockDim property. We need template for static shared memory allocation
             * @param a_rev a in  the reverse order
             * @param m size of a_rev
             * @param lookup precalc for compressed string b
             * @param vector_v compressed vector
             * @param carries of size m. Carry bits for each row of a_rev
             * @param n_small size of compressed string b (n_small = n / 32)
             * @param offset_row from bottom border
             * @param offset_col from left border
             */
            template<int Iter, int ChunkSize, int Width, int BlockDim>
            __global__ void
            hyyro_shared(int *a_rev, int m, Bit32 *lookup, Bit32 *vector_v, Bit32 *carries, int n_small, int offset_row,
                         int offset_col) {

                int warps_per_threadblock = BlockDim / 32;
                int per_warp_cols = (32 - 1) + Width;                             // per warp processed columns number
                int per_block_cols = Iter * per_warp_cols * warps_per_threadblock;// per block processed columns number
                int lane_id = threadIdx.x % 32;                                   // thread id within warp
                int warp_id = threadIdx.x / 32;                                   // warp id of current thread

                // Each warp manage own memory space of size ChunkSize.
                // By varying ChunkSize we make  tradeoffs between occupancy and number of load operations to shared memory.
                // To eliminate bank conflicts  we fill store  in column major order
                __shared__ Bit32 precalced_part[ChunkSize * warps_per_threadblock][32];

                // Each warp store in shared memory  part of processed vector for one Iter
                __shared__ Bit32 vector_part[warps_per_threadblock * per_warp_cols];


                // global positions

                // thread position
                int global_row = offset_row + lane_id + 32 * Iter * (warp_id + blockIdx.x * warps_per_threadblock);
                //
                int global_col = offset_col + lane_id + warp_id * Iter * per_warp_cols + blockIdx.x * per_block_cols;

                int precalc_offset = warp_id * ChunkSize;
                int shared_offset_v = warp_id * per_warp_cols;

                // within each iteration each warp process 32 rows elements and  (32 - 1) + Width columns
#pragma unroll
                for (int k = 0; k < Iter; k++, global_row += 32, global_col += per_warp_cols) {

                    // load carry bits
                    Bit32 carry_bit = (global_row < m) ? carries[global_row] : 0;

                    // load symbol a aka key
                    //load key for each row, no need for shared memory since we can use shfl_sync
                    int private_key = (global_row < m) ? a_rev[global_row] : 0;

                    //  load vector: each warp load  per_warp_cols elements
#pragma unroll
                    for (int i = 0; i < (per_warp_cols + 31) / 32; i++) {
                        if (lane_id + i * 32 < per_warp_cols)
                            vector_part[lane_id + i * 32 + shared_offset_v] = vector_v[global_col + 32 * i];
                    }

#pragma unroll
                    for (int it = 0; it < Width / ChunkSize; it++) {

                        // load precalc from shared memory
                        // warp first load row elements for thread 0, then for thread 1 ... thread 31
#pragma unroll
                        for (int i = 0; i < 32; i++) {
                            // broadcast current private key  for load chunk of values for it
                            int key = __shfl_sync(0xffffffff, private_key, i, 32);
                            //Note: no need for o  since vectors > n_small will be zeroed
                            int pos = i + it * ChunkSize + global_col + key * n_small;
                            if (lane_id < ChunkSize && pos < n_small) precalced_part[lane_id + ChunkSize *
                                                                                               warp_id][i] = lookup[pos];
                        }

                        __syncwarp();


#pragma unroll
                        for (int i = 0; i < ChunkSize; i++) {
                            Bit32 v = vector_part[lane_id + i + it * ChunkSize + shared_offset_v];
                            Bit32 precalc = precalced_part[i + precalc_offset][lane_id];
                            Bit32 p = precalc & v;


                            v += p;
                            carry_bit = v < p;
                            v += carry_bit;
                            //note case when carry bit is 1, and p =0 and v = 2**n - 1
                            carry_bit |= (v < p);

                            v = (p ^ v) | v;

                            __syncwarp();

                            //save back v
                            vector_part[lane_id + i + it * ChunkSize + shared_offset_v] = v;
                            __syncwarp();
                        }
                    }



                    // save carry bit
                    if (global_row < m) carries[global_row] = carry_bit;

                    // save vector
#pragma unroll
                    for (int i = 0; i < (per_warp_cols + 31) / 32; i++) {
                        if (lane_id + i * 32 < per_warp_cols)
                            vector_v[global_col + 32 * i] = vector_part[lane_id + i * 32 + shared_offset_v];
                    }
                    __syncwarp();
                }
            }


        }


        namespace sklansky {


            /**
             *
             *
             * @param a sum  of  two vectors consisting of 32 32bit values
             * @param b  basically sum of two vectors plus one for each summator
             * @param
             * @return  Update values for A and B and returns the  upper carry bit  for a value A
             */


            /**
             * Evaluate sum of two 1024 bit.
             * Each 1024 bit vector represents as 32 32 bit values.
             * Let $u=u_{0}u_{1}..u{31}$ and $v=v_{0}v_{1}..v{31}$ be 1024 bit vectors that we need to sum.
             * Then $a_{i} = u_{i} + v_{i}$, $b_{i} = a{i} + 1$ with possible overflow
             * Note that we need not to use __syncwarp since we use sync shuffle
             * @tparam withIf
             * @param a    sum  of  two
             * @param b   is a + 1
             * @param p
             * @param lane_id
             * @return
             */
            template<bool withIf>
            inline __device__ Bit32 sklansky_sum_reduction(Bit32 &u, Bit32 &v, Bit32 &carry, Bit32 &p, Bit32 &lane_id) {

                Bit32 a = u + v;
                Bit32 carry_a = a > u;
                a += carry;
                carry_a |= (a > u);
                Bit32 b = a + 1;
                carry_b = carry_a | b;





//TODO
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
                    put_other_value &= active_mask;

                    if (withIf) {

                        if (put_other_value) {
                            a = b;
                            with_no_carry = with_carry;
                        }

                    } else {

                        a = (a & (put_other_value - 1)) | ((-put_other_value) & b);
                        with_no_carry = (~put_other_value & with_no_carry) |
                                        (put_other_value &
                                         with_carry); // if put_other_value == 1 then we put with_no_carry <- with_carry;

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
            template<bool withIf>
            inline __device__ void
            kawanami_core(int *a_rev, Bit32 *lookup, Bit32 &vector, int n_small, int lane_id, Bit32 &own_carry,
                          int global_id_y) {

                int loc = 0; // current position of packed carries for rows
                int bit_pos = 0; // current position of lower bit in loc
                Bit32 carry = 0; // either 1 or 0 for the lower bit of big vector number, for others always 0
                Bit32 carry_pack_to_save = 0; // packed carry bits to save, only upper carry bit is saved
                Bit32 processing_carry_pack = (lane_id == 0) ? own_carry : 0; // only  lower adder can have carry

#pragma unroll
                for (int i = 0; i < 1024; i++) {

                    if (global_id_y < 0) {
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

                    Bit32 sum_v = p + vector + carry;// sum+1>p
                    Bit32 sum_v_inc = sum_v + 1;

                    carry_pack_to_save |= (sklansky_sum_reduction<withIf>(sum_v, sum_v_inc, p, lane_id) << bit_pos);
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
            __global__ void
            hyrro_kawanami_kernel_without_shared(int *a_rev, int m, Bit32 *lookup, Bit32 *vector_v, Bit32 *carries,
                                                 int n_small, int offset_x, int offset_y) {

                int warp_id = threadIdx.x / 32;
                int lane_id = threadIdx.x % 32;
                int columns_per_block =
                        BlockDim * Iter;// processed by threads within one block since each thread process Iter columns

                // position relative to global
                int global_id_x =
                        offset_x + lane_id + 32 * Iter * warp_id + blockIdx.x * columns_per_block; // global row pos
                int global_id_y = offset_y + H * 1024 * (1 + warp_id + blockIdx.x * (BlockDim / 32)) -
                                  1;// global col pos could be negative to adjust with non divisible parts

                int global_carry_id = global_id_y / 32 - lane_id;

#pragma unroll
                for (int i = 0; i < Iter; i++) {
                    Bit32 vector = vector_v[global_id_x + 32 * i];

#pragma unroll
                    for (int rep = 0; rep < H; rep++) {
                        //load packed carries from global memory
                        unsigned int own_carry = (global_carry_id - 32 * rep >= 0) ? carries[global_carry_id - 32 * rep]
                                                                                   : 0;
                        kawanami_core<WithIf>(a_rev, lookup, lane_id, own_carry, global_id_y);

                        //save carries for the 1024 elements
                        if (global_carry_id - 32 * rep >= 0) carries[global_carry_id] = own_carry;
                    }

                    // save 1024 bit vector back to memory
                    vector_v[global_id_x + 32 * i] = vector;
                }
            }


        }

    }
}


#endif //CPU_KAWANAMI_H
