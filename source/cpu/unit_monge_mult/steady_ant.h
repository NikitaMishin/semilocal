//
// Created by nikita on 15.10.2020.
//

#ifndef CPU_STEADY_ANT_H
#define CPU_STEADY_ANT_H

#include <algorithm>
#include <vector>
#include <iostream>
#include <climits>
// less then any value importand
#define NOPOINT (-1)


class Permutation {
private:
    int *row_to_col;
    int *col_to_row;

public:
    const int col_size; // rows aka height
    const int row_size; // cols aka width


    inline void set_point(int row, int col) {
        row_to_col[row] = col;
        col_to_row[col] = row;
    }

    inline void unset_point(int row, int col) {
        if (row_to_col[col] == col) {
            row_to_col[col] == NOPOINT;
            col_to_row[row] == NOPOINT;
        }
    }

    inline void unset_all() {
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }

    bool is_equal_to(Permutation &other) {
        if (other.row_size != row_size || other.col_size != col_size) return false;
        for (int i = 0; i < row_size; ++i) {
            if (get_col_by_row(i) != other.get_col_by_row(i)) return false;
        }
        for (int i = 0; i < col_size; ++i) {
            if (get_row_by_col(i) != other.get_row_by_col(i)) return false;
        }


        return true;

    }

    inline int get_row_by_col(int col) { return col_to_row[col]; }

    inline int get_col_by_row(int row) { return row_to_col[row]; }


    void print(std::ostream &os) {
        for (int i = 0; i < row_size; ++i) {
            for (int j = 0; j < col_size; ++j) {
                if (get_col_by_row(i) == j) {
                    os << 1;
                } else {
                    os << 0;
                }

            }
            os << std::endl;
        }
    }

    Permutation(int row, int col) : col_size(col), row_size(row) {
        row_to_col = new int[row];
        col_to_row = new int[col];
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }

    ~Permutation() {
        delete[] row_to_col;
        delete[] col_to_row;
    };
};

class Matrix {
private:
    int *arr;

public:
    const int row_size;
    const int col_size;

    Matrix(int row_size, int col_size) : row_size(row_size), col_size(col_size) {
        arr = new int[row_size * col_size];
        for (int i = 0; i < row_size * col_size; ++i) arr[i] = 0;
    }

    inline int get_element_at(int row, int col) { return arr[row * col_size + col]; }

    inline void set_element_at(int row, int col, int value) { arr[row * col_size + col] = value; }


    void print(std::ostream &os) {

        for (int i = 0; i < row_size; ++i) {
            for (int j = 0; j < col_size; ++j) {
                os << get_element_at(i, j);
            }
            os << std::endl;
        }
    }

    ~Matrix() {
        delete[] arr;
    }
};


namespace dominance_sum_counting {
    //ok
    namespace top_left_arrow {

        inline int left_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);
            if (row_cap >= row && row_cap != NOPOINT) sum++;
            return sum;
        }

        inline int down_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap >= col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_row_by_col(row - 1);
            if (col_cap >= col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);
            if (row_cap >= row && row_cap != NOPOINT) sum--;
            return sum;
        }


    }
    //ok
    namespace bottom_right_arrow {

        inline int left_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);

            if (row_cap != NOPOINT && row_cap < row) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int up_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int right_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);
            if (row_cap < row && row_cap != NOPOINT) sum++;
            return sum;
        }

    }

    namespace top_right_arrow {


        inline int left_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);
            if (row_cap >= row && row_cap != NOPOINT) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, Permutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);

            if (row_cap >= row && row_cap != NOPOINT) sum++;
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
    std::pair<int *, Permutation *> get_vertical_slice(Permutation *p_i, int col_start_inclusive, int col_exclusive) {
        auto new_size = col_exclusive - col_start_inclusive;

        auto cur_row_to_prev_row_mapping = new int[new_size];
        auto succ_pi = new Permutation(new_size, new_size);

        auto ptr = 0;
        for (int row = 0; row < p_i->row_size; ++row) {
            auto col = p_i->get_col_by_row(row);
            if (col >= col_start_inclusive & col < col_exclusive) {
                cur_row_to_prev_row_mapping[ptr] = row;
                succ_pi->set_point(ptr, col - col_start_inclusive);
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
    std::pair<int *, Permutation *> get_horizontal_slice(Permutation *q_i, int row_start_inclusive, int row_exclusive) {
        auto new_size = row_exclusive - row_start_inclusive;

        auto cur_col_to_prev_col_mapping = new int[new_size];
        auto succ_pi = new Permutation(new_size, new_size);

        auto ptr = 0;
        for (int col = 0; col < q_i->col_size; ++col) {
            auto row = q_i->get_row_by_col(col);
            if (row >= row_start_inclusive & row < row_exclusive) {
                cur_col_to_prev_col_mapping[ptr] = col;
                succ_pi->set_point(row - row_start_inclusive, ptr);

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
    inline void
    inverse_mapping(Permutation *shrinked, const int *row_mapper, const int *col_mapper, Permutation *flattened) {
        //could be parallelized
        flattened->unset_all();

        for (int cur_col = 0; cur_col < shrinked->col_size; ++cur_col) {
            auto old_col = col_mapper[cur_col]; /// consecutive access
            auto cur_row = shrinked->get_row_by_col(cur_col); // consecutive access
            auto old_row = row_mapper[cur_row]; // non consecutive access
            flattened->set_point(old_row, old_col);
        }
    }

    //TODO if use precalced matrices then what about free and so on?
    Permutation *steady_ant(Permutation *p, Permutation *q) {
        auto n = p->col_size;

        if (n == 1) {
            //base case when common dimension is one
            auto row = p->get_row_by_col(0);
            auto col = q->get_col_by_row(0);
            auto matrix = new Permutation(p->row_size, q->col_size);

            if (row != NOPOINT && col != NOPOINT) matrix->set_point(row, col);
            return matrix;
        }

        int spliter = p->col_size / 2;
        auto p_lo_tuple = get_vertical_slice(p, 0, spliter);
        auto p_hi_tuple = get_vertical_slice(p, spliter, p->col_size);


        auto q_lo_tuple = get_horizontal_slice(q, 0, spliter);

        auto r_lo = p; // reuse since we no need to use p further
        auto product = steady_ant(p_lo_tuple.second, q_lo_tuple.second);
        inverse_mapping(product, p_lo_tuple.first, q_lo_tuple.first, r_lo);
//        delete product;

        auto q_hi_tuple = get_horizontal_slice(q, spliter, q->row_size);
        auto r_hi = q; // reuse since we no need to use p further
        product = steady_ant(p_hi_tuple.second, q_hi_tuple.second);
        inverse_mapping(product, p_hi_tuple.first, q_hi_tuple.first, r_hi);
//        delete product;


        //ant passage
        auto end_row = -1;
        auto end_col = n + 1;
        auto cur_row = n;
        auto cur_col = -1;

        auto rhi = 0;
        auto rlo = 0;

        std::vector<int> good_row_pos;
        std::vector<int> good_col_pos;
//        good_row_pos.reserve(n / 2);
//        good_col_pos.reserve(n / 2);


        bool is_went_right = false; // went up
        while (true) {
//        for (int i = 0; i < n + n + 2; ++i) {

            if (end_col == cur_col && end_row == cur_row) break;

            if (cur_row == 0) break;
            //TODO is new
            if (cur_col == n) break;
            //
            auto dominance_row = cur_row - 1;
            auto dominance_col = cur_col + 1;

            //prev step
            if (is_went_right) {
                rhi = dominance_sum_counting::bottom_right_arrow::right_move(dominance_row, dominance_col - 1, rhi,
                                                                             r_hi);
                rlo = dominance_sum_counting::top_left_arrow::right_move(dominance_row, dominance_col - 1, rlo, r_lo);
            } else {
                rhi = dominance_sum_counting::bottom_right_arrow::up_move(dominance_row + 1, dominance_col, rhi, r_hi);
                rlo = dominance_sum_counting::top_left_arrow::up_move(dominance_row + 1, dominance_col, rlo, r_lo);
            }

            if (rhi - rlo < 0) {
                is_went_right = true;
                cur_col++;
            } else if (rhi - rlo == 0) {
                is_went_right = false;
                cur_row--;
            } else {
                std::cout << "Impissble" << std::endl;
            }

            if (dominance_col > 0) {
                auto delta_above_left =
                        dominance_sum_counting::bottom_right_arrow::left_move(dominance_row, dominance_col, rhi, r_hi) -
                        dominance_sum_counting::top_left_arrow::left_move(dominance_row, dominance_col, rlo, r_lo);
                auto delta_below_right =
                        dominance_sum_counting::bottom_right_arrow::down_move(dominance_row, dominance_col, rhi, r_hi) -
                        dominance_sum_counting::top_left_arrow::down_move(dominance_row, dominance_col, rlo, r_lo);

                if (delta_above_left < 0 && delta_below_right > 0) {
                    good_row_pos.push_back(dominance_row);
                    good_col_pos.push_back(dominance_col - 1);

                }
            }


        }
        // end ant passage
        std::cout << "R_LO and R_HI" << std::endl;
        r_lo->print(std::cout);
        std::cout << std::endl;
        r_hi->print(std::cout);
        std::cout << std::endl;
        // merge r_lo to r_hi
        for (int i = 0; i < n; ++i) {
            auto col = r_lo->get_col_by_row(i);
            if (col == NOPOINT) continue;
            r_hi->set_point(i, col);
        }

        // add good points
        for (int i = 0; i < good_col_pos.size(); ++i) {
            std::cout << "GOODPTS:";
            auto col = good_col_pos[i];
            auto row = good_row_pos[i];
            std::cout << "(" << row << "," << col << "),";
            r_hi->set_point(row, col);
        }
        std::cout << std::endl;

        return r_hi;
    }
}

namespace distance_product {

    Permutation *get_permutation_matrix(int row_size, int col_size, int seed = 0) {
        auto m = new Permutation(row_size, col_size);

        auto is_used = new bool[col_size];
        for (int i = 0; i < col_size; ++i) is_used[i] = false;

        /* initialize random seed: */
        srand(seed);

        auto active_pts = col_size;

        for (int row = 0; row < row_size && active_pts > 0; ++row) {

            while (true) {
                auto col = abs(rand()) % col_size;
                if (!is_used[col]) {
                    m->set_point(row, col);
                    is_used[col] = true;
                    break;
                }
            }
            active_pts--;
        }

        delete[] is_used;
        return m;
    }

    bool top_right_summator(int cur_row, int cur_col, int row_bound, int col_bound) {
        return cur_row >= row_bound && cur_col < col_bound;
    }


    namespace naive {

        template<typename Lambda>
        Matrix *get_dominance_matrix(Permutation &m, Lambda &&func) {
            auto row_size = m.row_size + 1;
            auto col_size = m.col_size + 1;
            auto dominance_matrix = new Matrix(row_size, col_size);

            for (int row = 0; row < row_size; ++row) {
                for (int col = 0; col < col_size; ++col) {
                    for (int i = 0; i < m.row_size; ++i) {
                        auto row_pos_point = i;
                        auto col_pos_point = m.get_col_by_row(row_pos_point);
                        if (col_pos_point == NOPOINT) continue;
                        if (func(row_pos_point, col_pos_point, row, col) == true)
                            dominance_matrix->set_element_at(row, col, dominance_matrix->get_element_at(row, col) + 1);
                    }
                }
            }
            return dominance_matrix;
        }


        Permutation *mult_dist(Permutation *m, Permutation *n) {
            auto dominance_m = get_dominance_matrix(*m, top_right_summator);
            auto dominance_n = get_dominance_matrix(*n, top_right_summator);
            auto row_size = (m->row_size + 1);
            auto col_size = (n->col_size + 1);
            auto dominance_c = new Matrix(row_size, col_size);

            for (int i = 0; i < row_size; ++i) {
                for (int k = 0; k < col_size; ++k) {
                    auto tmp = INT_MAX;
                    for (int j = 0; j < (m->col_size + 1); ++j) {
                        tmp = std::min(dominance_m->get_element_at(i, j) + dominance_n->get_element_at(j, k), tmp);
                    }
                    dominance_c->set_element_at(i, k, tmp);
                }
            }


            auto c = new Permutation(row_size - 1, col_size - 1);
            for (int i = 0; i < m->row_size; ++i) {
                for (int j = 0; j < n->col_size; ++j) {
                    auto cross_diff =
                            (dominance_c->get_element_at(i, j + 1) + dominance_c->get_element_at(i + 1, j)) -
                            (dominance_c->get_element_at(i, j) + dominance_c->get_element_at(i + 1, j + 1));

                    if (cross_diff == 1) c->set_point(i, j);

                }
            }

            delete dominance_c;
            delete dominance_m;
            delete dominance_n;

            return c;

        }
    }
}

#endif //CPU_STEADY_ANT_H
