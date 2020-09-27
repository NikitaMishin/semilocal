#include <cooperative_groups.h>

#define BRAID_ONES TODO
using namespace cooperative_groups; // or...
using cooperative_groups::thread_group; // etc.


template<class Input>
__device__ inline void prefix_braid_init(Input *arr, Input & value, int & pos) {
    arr[pos] = value;
}


template<class Input>
__device__ inline void process_cube_withoutif(Input & symbol_a, Input & left_strand, Input & symbol_b, Input & top_strand,
                                    Input & l_active_mask, Input & r_active_mask, bool & use_with_mask) {

    Input left_cap, symbols, combing_condition, rev_combing_cond, top_strand_shifted;

    int rev_counter = (sizeof(Input) * 8 - 2);
    Input mask = Input(1);
    Input braid_ones = Input(BRAID_ONES);

    // upper half
    #pragma unroll
    for (int rev_counter = (sizeof(Input) * 8 - 2); rev_counter > 0; rev_counter -= 2) {
        left_cap = left_strand >> rev_counter;
        symbols = ~(((symbol_a >> rev_counter)) ^ symbol_b);
        symbols &= (symbols >> 1) & braid_ones;
        combing_condition = mask & (symbols | (((~(left_cap)) & top_strand)));

        if(use_with_mask) {
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

    if(use_with_mask) {
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

        if(use_with_mask) {
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



// total thds = a_size




