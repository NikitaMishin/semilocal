#include <cooperative_groups.h>

using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.



template<class Input>
__device__ inline void prefix_braid_fill_cube_with_if(
        int left_edge, int top_edge, Input symbol_a_pack, Input symbol_b_pack,
        Input braid_ones,
        Input *bitset_left_strand_map,
        Input *bitset_top_strand_map,
        Input l_active_mask,
        Input r_active_mask,
        bool should_use_constraint) {

    Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;

    // access to  global date
    Input left_strand_pack = bitset_left_strand_map[left_edge];
    Input top_strand_pack = bitset_top_strand_map[top_edge];


    int rev_counter = (sizeof(Input) * 8 - 2);
    Input mask = Input(1);
    Input mask_r = Input(1) << rev_counter;


#pragma  unroll
    for (int i = 0; i < sizeof(Input) * 8 / 2; i++, rev_counter -= 2) {
        // upper fill
        left_cap = left_strand_pack >> rev_counter;
        symbols = ~(((symbol_a_pack >> rev_counter)) ^ symbol_b_pack);
        symbols &= (symbols >> 1) & braid_ones;

        combing_condition = mask & (symbols | (((~(left_cap)) & top_strand_pack)));

        if (should_use_constraint) {
            combing_condition &= (l_active_mask >> rev_counter) & r_active_mask;
        }

        rev_combing_cond = combing_condition ^ braid_ones;

        // its registers!
        if (combing_condition) {
            top_strand_shifted = top_strand_pack << rev_counter;
            top_strand_pack = (rev_combing_cond & top_strand_pack) | (combing_condition & left_cap);

            symbols = ~(((symbol_a_pack)) ^ (symbol_b_pack << rev_counter));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = mask_r & (symbols | ((~(left_strand_pack) & top_strand_shifted)));

            if (should_use_constraint) {
                combing_condition &= l_active_mask & (r_active_mask << rev_counter);
            }

            rev_combing_cond = combing_condition ^ braid_ones;

            left_strand_pack = (rev_combing_cond & left_strand_pack) | (combing_condition & top_strand_shifted);
        }

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);


    }

    // middle
    symbols = (~(symbol_a_pack ^ symbol_b_pack));
    symbols &= (symbols >> 1) & braid_ones;
    combing_condition = (symbols | ((~left_strand_pack) & top_strand_pack));

    if (should_use_constraint) {
        combing_condition &= (l_active_mask & r_active_mask);
    }

    rev_combing_cond = combing_condition ^ braid_ones;

    if (combing_condition) {
        top_strand_shifted = top_strand_pack;
        top_strand_pack =
                (rev_combing_cond & top_strand_pack) | (combing_condition & left_strand_pack);
        left_strand_pack =
                (rev_combing_cond & left_strand_pack) | (combing_condition & top_strand_shifted);
    }


    mask = braid_ones;
    mask_r = braid_ones;

    // lower
#pragma  unroll
    for (int i = 0; i < sizeof(Input) * 8 / 2; i++) {
        mask <<= 2;
        mask_r >>= 2;

        left_cap = left_strand_pack << (2 * (i + 1));
        symbols = ~(((symbol_a_pack << (2 * (i + 1)))) ^ symbol_b_pack);
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = mask & (symbols | (((~(left_cap)) & top_strand_pack)));

        if (should_use_constraint) {
            combing_condition &= (l_active_mask << (2 * (i + 1))) & r_active_mask;
        }

        rev_combing_cond = combing_condition ^ braid_ones;

        if (combing_condition) {
            top_strand_shifted = top_strand_pack >> (2 * (i + 1));

            top_strand_pack = (rev_combing_cond & top_strand_pack) | (combing_condition & left_cap);

            symbols = ~(((symbol_a_pack)) ^ (symbol_b_pack >> (2 * (i + 1))));
            symbols &= (symbols >> 1) & braid_ones;
            combing_condition = mask_r & (symbols | ((~(left_strand_pack) & top_strand_shifted)));

            if (should_use_constraint) {
                combing_condition &= l_active_mask & (r_active_mask >> (2 * (i + 1)));
            }

            rev_combing_cond = combing_condition ^ braid_ones;

            left_strand_pack = (rev_combing_cond & left_strand_pack) | (combing_condition & top_strand_shifted);
        }
    }

    // assign to global
    bitset_left_strand_map[left_edge] = left_strand_pack;
    bitset_top_strand_map[top_edge] = top_strand_pack;

}


template<class Input>
__device__ inline void prefix_braid_fill_cube_without_if(
        int left_edge, int top_edge, Input symbol_a_pack, Input symbol_b_pack, Input braid_ones,
        Input *bitset_left_strand_map,
        Input *bitset_top_strand_map,
        Input l_active_mask,
        Input r_active_mask,
        bool should_use_constraint) {

    Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;

    // access to  global date
    Input left_strand_pack = bitset_left_strand_map[left_edge];
    Input top_strand_pack = bitset_top_strand_map[top_edge];


    int rev_counter = (sizeof(Input) * 8 - 2);
    Input mask = Input(1);
    Input mask_r = Input(1) << rev_counter;


#pragma  unroll
    for (int i = 0; i < sizeof(Input) * 8 / 2; i++, rev_counter -= 2) {
        // upper fill
        left_cap = left_strand_pack >> rev_counter;
        symbols = ~(((symbol_a_pack >> rev_counter)) ^ symbol_b_pack);
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = mask & (symbols | (((~(left_cap)) & top_strand_pack)));

        if (should_use_constraint) {
            combing_condition &= (l_active_mask >> rev_counter) & r_active_mask;
        }

        rev_combing_cond = combing_condition ^ braid_ones;

        // its registers!
        top_strand_shifted = top_strand_pack << rev_counter;
        top_strand_pack = (rev_combing_cond & top_strand_pack) | (combing_condition & left_cap);

        symbols = ~(((symbol_a_pack)) ^ (symbol_b_pack << rev_counter));
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = mask_r & (symbols | ((~(left_strand_pack) & top_strand_shifted)));

        if (should_use_constraint) {
            combing_condition &= l_active_mask & (r_active_mask << rev_counter);
        }


        rev_combing_cond = combing_condition ^ braid_ones;

        left_strand_pack = (rev_combing_cond & left_strand_pack) | (combing_condition & top_strand_shifted);

        mask = (mask << 2) | Input(1);
        mask_r = mask_r | (mask_r >> 2);


    }

    // middle
    symbols = (~(symbol_a_pack ^ symbol_b_pack));
    symbols &= (symbols >> 1) & braid_ones;
    combing_condition = (symbols | ((~left_strand_pack) & top_strand_pack));

    if (should_use_constraint) {
        combing_condition &= (l_active_mask & r_active_mask);
    }


    rev_combing_cond = combing_condition ^ braid_ones;

    top_strand_shifted = top_strand_pack;

    top_strand_pack =
            (rev_combing_cond & top_strand_pack) | (combing_condition & left_strand_pack);
    left_strand_pack =
            (rev_combing_cond & left_strand_pack) | (combing_condition & top_strand_shifted);


    mask = braid_ones;
    mask_r = braid_ones;

    // lower
#pragma  unroll
    for (int i = 0; i < sizeof(Input) * 8 / 2; i++) {
        mask <<= 2;
        mask_r >>= 2;

        left_cap = left_strand_pack << (2 * (i + 1));
        symbols = ~(((symbol_a_pack << (2 * (i + 1)))) ^ symbol_b_pack);
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = mask & (symbols | (((~(left_cap)) & top_strand_pack)));

        if (should_use_constraint) {
            combing_condition &= (l_active_mask << (2 * (i + 1))) & r_active_mask;
        }


        rev_combing_cond = combing_condition ^ braid_ones;

        top_strand_shifted = top_strand_pack >> (2 * (i + 1));

        top_strand_pack = (rev_combing_cond & top_strand_pack) | (combing_condition & left_cap);

        symbols = ~(((symbol_a_pack)) ^ (symbol_b_pack >> (2 * (i + 1))));
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = mask_r & (symbols | ((~(left_strand_pack) & top_strand_shifted)));

        if (should_use_constraint) {
            combing_condition &= l_active_mask & (r_active_mask >> (2 * (i + 1)));
        }

        rev_combing_cond = combing_condition ^ braid_ones;

        left_strand_pack = (rev_combing_cond & left_strand_pack) | (combing_condition & top_strand_shifted);

    }

    // assign to global
    bitset_left_strand_map[left_edge] = left_strand_pack;
    bitset_top_strand_map[top_edge] = top_strand_pack;
}

__device__ inline void
semi_local_init_phase_conseq_fill_block(int *reduced_sticky_braid, int thread_id_global, int offset) {
    reduced_sticky_braid[thread_id_global] = offset + thread_id_global;
}


template<class Input>
__device__ inline void prefix_braid_init(Input *arr, int arr_size, Input value, int pos) {
    arr[pos] = value;
}


// total thds = a_size
template<class Input>
__global__ void semi_local_gpu_withif(Input const *seq_a_rev, int a_size, Input const *seq_b, int b_size,
                                      Input *bitset_left_strand_map, Input *bitset_top_strand_map,
                                      int cells_per_thread, int total_thds, Input braid_ones, Input l_active,
                                      Input r_active) {

    // a < b
    //thread_id = m

    // is this thread
    Input l_active_pack, r_active_pack;
    bool should = false;


    int num_diag = a_size + b_size - 1;
    int total_same_length_diag = num_diag - (a_size - 1) - (a_size - 1);
    // in [0..a_size-1]
    auto thread_id = blockIdx.x * blockDim.x + threadIdx.x;


    //sync primitive to sync whole grid
    auto g = this_grid();

    //init_phase
    for (int i = 0; i < cells_per_thread; i++) {
        if (total_thds * i + thread_id < a_size + b_size)
            prefix_braid_init(bitset_left_strand_map, braid_ones, thread_id + total_thds * i);
    }



    g.sync();

    // phase one and three
    for (int i = 0; i < a_size - 1; i++) {
        //todo check and inline
        int is_first_phase = thread_id >= a_size - 1 - 1 - i;
        if () {
            // thread_id
            // phase one
            //modulo
            prefix_braid_fill_cube_without_if();
        } else {

            prefix_braid_fill_cube_without_if();
        }
        g.sync();
    }


    //phase two
    for (int i = 0; i < total_same_length_diag; i++) {
        prefix_braid_fill_cube_without_if();
        g.sync();
    }

}





