

#include <cooperative_groups.h>


namespace semi_local_lcs_gpu {
    using namespace cooperative_groups;
    using cooperative_groups::thread_group;

    namespace gpu_routines {


        /***
        * Initial phase. Filling braid:
        * a[i] = i \forall i \in [0...a.size())
        * @param reduced_sticky_braid
        * @param thread_id_global
        * @param offset
        */
        template<class Input>
        __device__ inline void
        semi_local_init_phase_conseq_fill(Input *reduced_sticky_braid, Input pos) {
            reduced_sticky_braid[pos] = pos;
        }

        /**
         * Process diagonal
         * @tparam Input
         * @param l_strand
         * @param top_strand
         * @param a_symbol
         * @param seq_b
         * @param top
         * @param b_pos
         */
        template<class Input>
        __device__ inline void
        semi_local_diag_fill(Input &l_strand, Input &top_strand, Input &a_symbol, Input &b_symbol) {

            int old_left = l_strand;
            int should_swap = (a_symbol == b_symbol) || (l_strand > top_strand);

            // aka swap if true else id
            l_strand = (1 - should_swap) * l_strand + should_swap * top_strand;
            top_strand = (1 - should_swap) * top_strand + should_swap * old_left;
        }


        /**
         *
         * @tparam Input type holder of elements of sticky braid and sequence a and b
         * @param sticky_braid
         * @param seq_a
         * @param a_size
         * @param seq_b
         * @param b_size
         * @param offset_l
         * @param offset_top
         * @param diag_len
         * @param cells_per_thd
         * @param total_thds
         */
        template<class Input>
        __global__ void
        process_diagonal(Input *sticky_braid, Input const *seq_a, int a_size, Input const *seq_b, int b_size,
                         int offset_l,
                         int offset_top, int diag_len, int cells_per_thd, int total_thds) {

            for (int i = 0; i < cells_per_thd; i++) {
                int thread_id = blockIdx.x * blockDim.x + threadIdx.x + i * total_thds;

                if (thread_id < diag_len) {
                    //load left symbol
                    Input a_symbol = seq_a[a_size - offset_l - 1 - thread_id];
                    Input b_symbol = seq_b[thread_id + offset_top];
                    Input left_strand = sticky_braid[thread_id + offset_l];
                    Input top_strand = sticky_braid[thread_id + offset_top + a_size];
                    semi_local_diag_fill(left_strand, top_strand, a_symbol, b_symbol);

                    sticky_braid[thread_id + offset_l] = left_strand;
                    sticky_braid[thread_id + offset_top + a_size] = top_strand;
                }
            }


        }



//        __global__ void semi_local_withoutif_prestored_lefts(int const *seq_a, int a_size, int const *seq_b, int b_size,
//                                                             int *reduced_sticky_braid, int cells_per_thread, int total_thds) {
//
//
//            int num_diag = a_size + b_size - 1;
//            int total_same_length_diag = num_diag - (a_size - 1) - (a_size - 1);
//            int thread_id = blockIdx.x * blockDim.x + threadIdx.x;
//
//            //sync primitive to sync whole grid
//            auto g = this_grid();
//
//            //init_phase
//            for (int i = 0; i < cells_per_thread; i++) {
//                if (total_thds * i + thread_id < a_size + b_size)
//                    semi_local_init_phase_conseq_fill(reduced_sticky_braid, thread_id + total_thds * i);
//            }
//
//
//            // todo do we need apply forward access pattern and then swap within warp?
//            int left_symbol = (a_size - 1 - thread_id) >= 0 ? seq_a[a_size - 1 - thread_id] : 0;
//            // economy of
//            int left_strand = thread_id;
//
//            g.sync();
//
//            //1 phase
//            int b_pos = 0;
//            for (int i = a_size - 1; i > 0; i--) {
//                // only specific threads  && only active thread should perform
//                if (thread_id >= i && thread_id < a_size) {
//                    int top_strand = reduced_sticky_braid[b_pos + a_size];
//                    semi_local_diag_fill(reduced_sticky_braid, left_strand,top_strand, left_symbol, seq_b, a_size, b_pos);
//                    reduced_sticky_braid[b_pos + a_size] = top_strand;
//                    b_pos++;
//                }
//                g.sync();
//            }
//
//            //2 phase
//            for (int i = 0; i < total_same_length_diag; i++) {
//                if (thread_id < a_size) {
//                    int top_strand = reduced_sticky_braid[b_pos + a_size];
//                    semi_local_diag_fill(reduced_sticky_braid, left_strand,top_strand, left_symbol, seq_b, a_size, b_pos);
//                    reduced_sticky_braid[b_pos + a_size] = top_strand;
//                    b_pos++;
//                }
//                g.sync();
//            }
//
//            //3 phase
//            for (int i = a_size - 2; i >= 0; i--) {
//                if (thread_id <= i) {
//                    int top_strand = reduced_sticky_braid[b_pos + a_size];
//                    semi_local_diag_fill(reduced_sticky_braid, left_strand,top_strand, left_symbol, seq_b, a_size, b_pos);
//                    reduced_sticky_braid[b_pos + a_size] = top_strand;
//                    b_pos++;
//                }
//                g.sync();
//            }
//
//            reduced_sticky_braid[thread_id] = left_strand;
//        }



    }


    namespace gpu_wrappers {

        /**
         * TODO
         * @tparam Input
         * @param a_gpu
         * @param a_size
         * @param b_gpu
         * @param b_size
         * @param sticky_braid_gpu
         * @param block_size
         * @param total_thds
         * @param offset_l
         * @param offset_top
         * @param diag_len
         */
        template<class Input>
        void semi_local_process_diagonal(Input *a_gpu, int a_size, Input *b_gpu, int b_size,
                                         Input *sticky_braid_gpu, int block_size, int total_thds,
                                         int offset_l, int offset_top, int diag_len) {


            int total_blocks = std::ceil(float(total_thds) / block_size);
            total_thds = total_blocks * block_size;
            int cells_per_thd = std::ceil(float(diag_len) / total_thds);

            dim3 grid(total_blocks, 1);
            dim3 block(block_size, 1);


            std::cout << cells_per_thd << std::endl;
            std::cout << total_thds << std::endl;
            std::cout << a_size << std::endl;
            std::cout << b_size << std::endl;
            std::cout << diag_len << std::endl;
            std::cout << "s";


            gpu_routines::process_diagonal<<<grid, block>>>(sticky_braid_gpu, a_gpu, a_size, b_gpu, b_size,
                    offset_l, offset_top, diag_len, cells_per_thd, total_thds);
            memory_management::gpuAssert(cudaGetLastError(), __FILE__, __LINE__);
            memory_management::synchronize_with_gpu();

        }

        void semi_local_lcs_process_all_diagonals(){
//            TODO();
        }



    }

}












