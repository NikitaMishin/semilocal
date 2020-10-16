//
// Created by nikita on 15.10.2020.
//

#ifndef CPU_STEADY_ANT_H
#define CPU_STEADY_ANT_H

#include <algorithm>
#include <vector>
#include <iostream>
// less then any value importand
#define NOPOINT (-1)

struct Matrix {
    int *row_to_col;
    int *col_to_row;
    int col_size; // rows aka height
    int row_size; // cols aka width
};

namespace dominance_sum_counting {

    namespace top_left_arrow {

        inline int left_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->col_to_row[col - 1];
            if (row_cap >= row) sum++;
            return sum;
        }

        inline int down_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->row_to_col[row];
            if (col_cap >= col) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->row_to_col[row - 1];
            if (col_cap >= col) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->col_to_row[col];
            if (row_cap >= row) sum--;
            return sum;
        }


    }

    namespace bottom_right_arrow {

        inline int left_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->col_to_row[col - 1];

            if (row_cap != NOPOINT && row_cap < row) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->row_to_col[row];
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int up_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->row_to_col[row - 1];
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int right_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->col_to_row[col];
            if (row_cap < row && row_cap != NOPOINT) sum++;
            return sum;
        }

    }

    namespace top_right_arrow {


        inline int left_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->col_to_row[col - 1];
            if (row_cap >= row) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->row_to_col[row];
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->row_to_col[row - 1];
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, Matrix *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->col_to_row[col];

            if (row_cap >= row) sum++;
            return sum;
        }

    }

};

namespace steady_ant {



    //TODO maybe use sorting algo for slicing? we now range u..k of keys, we need just sort associated values to place in 0..k-u
    // TODO really case when parallel will be


    /**
     * Given squared permutation matrix p_{i} get slice p[,start_inclusive:end_exclusive] and map it to new coordinates
     * to get new squared matrix of size end_exclusive-start_inclusive and mapping of row coordinates (new to old)
     * @param p_i squared permutation matrix
     * @param col_start_inclusive inclusive index of start of slice
     * @param col_exclusive exclusive index of end of slice
     * @return mapping of row_{i+1} -> row_{i} and p_{i+1}, mapping of col_{i+1} -> col_{i} implicit via offset start_inclusive
     */
    std::pair<int *, Matrix *> get_vertical_slice(Matrix *p_i, int col_start_inclusive, int col_exclusive) {
        auto new_size = col_exclusive - col_start_inclusive;

        auto cur_row_to_prev_row_mapping = new int[new_size];
        auto succ_pi = new Matrix();
        succ_pi->row_size = new_size;
        succ_pi->col_size = new_size;
        succ_pi->row_to_col = new int[new_size];
        succ_pi->col_to_row = new int[new_size];

        auto ptr = 0;
        for (int row = 0; row < p_i->row_size; ++row) {
            auto col = p_i->row_to_col[row];
            if (col >= col_start_inclusive & col < col_exclusive) {
                cur_row_to_prev_row_mapping[ptr] = row;
                succ_pi->row_to_col[ptr] = col - col_start_inclusive;
                succ_pi->col_to_row[col - col_start_inclusive] = ptr;

                ptr++;
            }
        }

        return std::make_pair(cur_row_to_prev_row_mapping, succ_pi);
    }


    /**
     * Given squared permutation matrix q_{i} get slice p[start_inclusive:end_exclusive,:] and map it to new coordinates
    * to get new squared matrix of size end_exclusive-start_inclusive and mapping of col coordinates (new to old)
    * @param q_i squared permutation matrix
    * @param row_start_inclusive inclusive index of start of slice
    * @param row_exclusive exclusive index of end of slice
    * @return mapping of col_{i+1} -> col_{i} and p_{i+1}, mapping of row_{i+1} -> row_{i} implicit via offset start_inclusive
    */
    std::pair<int *, Matrix *> get_horizontal_slice(Matrix *q_i, int row_start_inclusive, int row_exclusive) {
        auto new_size = row_exclusive - row_start_inclusive;

        auto cur_col_to_prev_col_mapping = new int[new_size];
        auto succ_pi = new Matrix();
        succ_pi->row_size = new_size;
        succ_pi->col_size = new_size;
        succ_pi->row_to_col = new int[new_size];
        succ_pi->col_to_row = new int[new_size];

        auto ptr = 0;
        for (int col = 0; col < q_i->col_size; ++col) {
            auto row = q_i->col_to_row[col];
            if (row >= row_start_inclusive & col < row_exclusive) {
                cur_col_to_prev_col_mapping[ptr] = col;
                succ_pi->row_to_col[ptr] = row - row_start_inclusive;
                succ_pi->col_to_row[row - row_start_inclusive] = ptr;

                ptr++;
            }
        }

        return std::make_pair(cur_col_to_prev_col_mapping, succ_pi);
    }

    /**
     * also set -1 to no point
     * @param shrinked
     * @param row_mapper
     * @param col_mapper
     * @param flattened
     */
    inline void inverse_mapping(Matrix *shrinked, const int *row_mapper, const int *col_mapper, Matrix *flattened) {
        //could be parallelized
        for (int i = 0; i < flattened->col_size; ++i) {
            flattened->row_to_col[i] = -1;
            flattened->col_to_row[i] = -1;

        }

        for (int cur_col = 0; cur_col < shrinked->col_size; ++cur_col) {
            auto old_col = col_mapper[cur_col]; /// consecutive access
            auto cur_row = shrinked->col_to_row[cur_col]; // consecutive access
            auto old_row = row_mapper[cur_row]; // non consecutive access
            flattened->col_to_row[old_col] = old_row;
            flattened->row_to_col[old_row] = old_col;
        }
    }

    //TODO if use precalced matrices then what about free and so on?
    Matrix *steady_ant(Matrix *p, Matrix *q) {
        auto n = p->col_size;

        if (n == 1) {
            //base case when common dimension is one
            auto row = p->col_to_row[0];
            auto col = q->row_to_col[0];
            auto matrix = new Matrix();
            matrix->row_size = p->row_size;
            matrix->col_size = q->col_size;
            matrix->row_to_col = new int[matrix->row_size];
            matrix->col_to_row = new int[matrix->col_size];

            for (int i = 0; i < matrix->row_size; ++i) {
                matrix->row_to_col[i] = NOPOINT;
            }
            for (int i = 0; i < matrix->col_size; ++i) {
                matrix->col_to_row[i] = NOPOINT;
            }

            if (row != NOPOINT && col != NOPOINT) {
                matrix->col_to_row[col] = row;
                matrix->row_to_col[row] = col;
            }

            return matrix;
        }

        int spliter = p->col_size / 2;
        auto p_lo_tuple = get_vertical_slice(p, 0, spliter);
        auto p_hi_tuple = get_vertical_slice(p, spliter, p->col_size);
        auto q_lo_tuple = get_horizontal_slice(q, 0, spliter);

        auto r_lo = p; // reuse since we no need to use p further
        auto product = steady_ant(p_lo_tuple.second, q_lo_tuple.second);
        inverse_mapping(product, p_lo_tuple.first, q_lo_tuple.first, r_lo);
        free(product);

        auto q_hi_tuple = get_horizontal_slice(q, spliter, q->row_size);//TODO row_size?
        auto r_hi = q; // reuse since we no need to use p further
        product = steady_ant(p_hi_tuple.second, q_hi_tuple.second);
        inverse_mapping(product, p_hi_tuple.first, q_hi_tuple.first, r_hi);
        free(product);


        //ant passage
        auto end_row = -1;
        auto end_col = n + 1;
        auto cur_row = r_lo->row_size;
        auto cur_col = -1;

        auto rhi = 0;
        auto rlo = 0;

        std::vector<int> good_row_pos;
        std::vector<int> good_col_pos;
        auto good_points = 0;
        good_row_pos.reserve(n / 2);
        good_col_pos.reserve(n / 2);


        bool is_went_right = false; // went up
        for (int i = 0; i < n + n + 2; ++i) {

            if (end_col == cur_col && end_row == cur_row) break;

            if (cur_row == 0) break;
            //TODO is new
            if (cur_col == n) break;
            //
            auto dominance_row = cur_row - 1;
            auto dominance_col = cur_col + 1;

            if (is_went_right) {
                rhi = dominance_sum_counting::bottom_right_arrow::right_move(dominance_row, cur_col, rhi, r_hi);
                rlo = dominance_sum_counting::top_left_arrow::right_move(dominance_row, cur_col, rlo, r_lo);
            } else {
                rhi = dominance_sum_counting::bottom_right_arrow::up_move(cur_row, dominance_col, rhi, r_hi);
                rlo = dominance_sum_counting::top_left_arrow::up_move(cur_row, dominance_col, rlo, r_lo);
            }

            if (r_hi - r_lo < 0) {
                is_went_right = true;
                cur_col++;
            } else if (r_hi - r_lo == 0) {
                is_went_right = false;
                cur_row--;
            } else {
                std::cout << "Impissble" << std::endl;
            }

            if (cur_col + 1 > 0) {
                auto delta_above_left =
                        dominance_sum_counting::bottom_right_arrow::left_move(dominance_row, dominance_col, rhi, r_hi) -
                        dominance_sum_counting::top_left_arrow::left_move(dominance_row, dominance_col, rlo, r_lo);
                auto delta_below_right =
                        dominance_sum_counting::bottom_right_arrow::down_move(dominance_row, dominance_col, rhi, r_lo) -
                        dominance_sum_counting::top_left_arrow::down_move(dominance_row, dominance_col, rlo, r_lo);

                if (delta_above_left < 0 && delta_below_right > 0) {
                    good_row_pos[good_points] = dominance_row;
                    good_col_pos[good_points] = dominance_col - 1;
                    good_points++;
                }
            }


        }
        // end ant passage

        // merge l_lo to r_hi
        for (int i = 0; i < n; ++i) {
            auto col = r_lo->row_to_col[i];
            if (col == -1) continue;
            r_hi->row_to_col[i] = col;
            r_hi->col_to_row[col] = i;
        }

        // add good points
        for (int i = 0; i < good_points; ++i) {
            auto col = good_col_pos[i];
            auto row = good_row_pos[i];
            r_hi->row_to_col[row] = col;
            r_hi->col_to_row[col] = row;
        }

        return r_hi;
    }
}


#endif //CPU_STEADY_ANT_H
