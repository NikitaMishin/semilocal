


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



