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

typedef std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>> PrecalcMap;


class AbstractPermutation {
public:
    int col_size; // rows aka height
    int row_size; // cols aka width

    inline virtual void set_point(int row, int col) = 0;

    inline virtual void unset_point(int row, int col) = 0;

    inline virtual void unset_all() = 0;


    inline virtual int get_row_by_col(int col) const = 0;

    inline virtual int get_col_by_row(int row) const = 0;

    virtual ~AbstractPermutation() = default;;


    AbstractPermutation(int row, int col) : col_size(col), row_size(row) {}


    void print(std::ostream &os) const {
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


    bool is_equal_to(const AbstractPermutation &other) const {
        if (other.row_size != row_size || other.col_size != col_size) return false;
        for (int i = 0; i < row_size; ++i) {
            if (get_col_by_row(i) != other.get_col_by_row(i)) return false;
        }
        for (int i = 0; i < col_size; ++i) {
            if (get_row_by_col(i) != other.get_row_by_col(i)) return false;
        }
        return true;
    }


    void to_points_on_grid(std::vector<std::pair<int, int>> &result) const {
        for (int i = 0; i < col_size; ++i) {
            auto col = get_col_by_row(i);
            if (col != NOPOINT) result.emplace_back(i, col);
        }
    }


};


class Permutation : public AbstractPermutation {
private:
    int *row_to_col;
    int *col_to_row;

public:

    inline void set_point(int row, int col) override {
        row_to_col[row] = col;
        col_to_row[col] = row;
    }

    inline void unset_point(int row, int col) override {
        if (row_to_col[col] == col) {
            row_to_col[col] = NOPOINT;
            col_to_row[row] = NOPOINT;
        }
    }

    inline void unset_all() override {
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }

    inline int get_row_by_col(int col) const override { return col_to_row[col]; }

    inline int get_col_by_row(int row) const override { return row_to_col[row]; }


    Permutation(int row, int col) : AbstractPermutation(row, col) {
        row_to_col = new int[row];
        col_to_row = new int[col];
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }

    Permutation(int row, int col, std::vector<std::pair<int, int>> &points) : Permutation(row, col) {
        for (auto &point: points) {
            row_to_col[point.first] = point.second;
            col_to_row[point.second] = point.first;
        }
    }

    ~Permutation() override {
        delete[] row_to_col;
        delete[] col_to_row;
    };


};

class PermutationPreAllocated : public AbstractPermutation {

public:

    int *row_to_col;
    int *col_to_row;

    inline void set_point(int row, int col) override {
        row_to_col[row] = col;
        col_to_row[col] = row;
    }

    inline void unset_point(int row, int col) override {
        if (row_to_col[col] == col) {
            row_to_col[col] = NOPOINT;
            col_to_row[row] = NOPOINT;
        }
    }

    inline void unset_all() override {
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
    }


    inline int get_row_by_col(int col) const { return col_to_row[col]; }

    inline int get_col_by_row(int row) const { return row_to_col[row]; }


    PermutationPreAllocated(int row, int col, int *row_to_col, int *col_to_row) : AbstractPermutation(row, col),
                                                                                  row_to_col(row_to_col),
                                                                                  col_to_row(col_to_row) {};


    PermutationPreAllocated(int row, int col, int *row_to_col, int *col_to_row,
                            std::vector<std::pair<int, int>> &points) :
            PermutationPreAllocated(row, col, row_to_col, col_to_row) {
        for (int i = 0; i < row_size; ++i) row_to_col[i] = NOPOINT;
        for (int i = 0; i < col_size; ++i) col_to_row[i] = NOPOINT;
        for (auto &point: points) {
            row_to_col[point.first] = point.second;
            col_to_row[point.second] = point.first;
        }
    }

    ~PermutationPreAllocated() override = default;;

};


namespace std {

    template<>
    struct hash<std::pair<int, int>> {
        std::size_t operator()(const std::pair<int, int> &pair) const {
            return hash<long long>()(((long long) pair.first) ^ (((long long) pair.second) << 32));
        }
    };


    template<>
    struct hash<AbstractPermutation> {
        std::size_t operator()(const AbstractPermutation &k) const {
            using std::size_t;
            using std::hash;
            using std::string;

            auto pass_by_row = k.row_size < k.col_size;
            size_t sum = 0;
            auto bits_per_symbol = int(std::ceil(log2(k.row_size)));

            for (int i = 0; i < k.row_size; ++i) {
                sum = (sum << bits_per_symbol) + k.get_col_by_row(i);

            }
            return sum;
        }

    };

    bool operator==(const Permutation &p1, const Permutation &p2) {
        return std::hash<AbstractPermutation>()(p1) == std::hash<AbstractPermutation>()(p1);
    }

}

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

        inline int left_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);
            if (row_cap >= row && row_cap != NOPOINT) sum++;
            return sum;
        }

        inline int down_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap >= col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap >= col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);
            if (row_cap >= row && row_cap != NOPOINT) sum--;
            return sum;
        }


    }
    //ok
    namespace bottom_right_arrow {

        inline int left_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);

            if (row_cap != NOPOINT && row_cap < row) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int up_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int right_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);
            if (row_cap < row && row_cap != NOPOINT) sum++;
            return sum;
        }

    }

    namespace top_right_arrow {


        inline int left_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col == 0) return sum;
            auto row_cap = perm_matrix->get_row_by_col(col - 1);
            if (row_cap >= row && row_cap != NOPOINT) sum--;
            return sum;
        }

        inline int down_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row >= perm_matrix->row_size) return 0;
            auto col_cap = perm_matrix->get_col_by_row(row);
            if (col_cap < col && col_cap != NOPOINT) sum--;
            return sum;
        }

        inline int up_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (row == 0) return sum;
            auto col_cap = perm_matrix->get_col_by_row(row - 1);
            if (col_cap < col && col_cap != NOPOINT) sum++;
            return sum;
        }

        inline int right_move(int row, int col, int sum, AbstractPermutation *perm_matrix) {
            if (col >= perm_matrix->col_size) return 0;
            auto row_cap = perm_matrix->get_row_by_col(col);

            if (row_cap >= row && row_cap != NOPOINT) sum++;
            return sum;
        }

    }

};

namespace distance_product {

    AbstractPermutation *get_permutation_matrix(int row_size, int col_size, int seed = 0) {
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
        Matrix *get_dominance_matrix(AbstractPermutation &m, Lambda &&func) {
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


        Permutation *mult_dist(AbstractPermutation *m, AbstractPermutation *n) {
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


namespace steady_ant {

    /**
     * Given squared permutation matrix p_{i} get slice p[,start_inclusive:end_exclusive] and map it to new coordinates
     * to get new squared matrix of size end_exclusive-start_inclusive and mapping of row coordinates (new to old)
     * @param p_i squared permutation matrix
     * @param col_start_inclusive inclusive index of start of slice
     * @param col_end_exclusive exclusive index of end of slice
     * @return mapping of row_{i+1} -> row_{i} and p_{i+1}, mapping of col_{i+1} -> col_{i} implicit via offset start_inclusive
     */
    void get_vertical_slice(AbstractPermutation *p_i, int col_start_inclusive, int col_end_exclusive,
                            int *cur_row_to_prev_row_mapping, AbstractPermutation *succ_pi) {

        auto ptr = 0;
        for (int row = 0; row < p_i->row_size; ++row) {
            auto col = p_i->get_col_by_row(row);
            if (col >= col_start_inclusive & col < col_end_exclusive) {
                cur_row_to_prev_row_mapping[ptr] = row;
                succ_pi->set_point(ptr, col - col_start_inclusive);
                ptr++;
            }
        }
    }


    /**
     * Given squared permutation matrix q_{i} get slice p[start_inclusive:end_exclusive,:] and map it to new coordinates
    * to get new squared matrix of size end_exclusive-start_inclusive and mapping of col coordinates (new to old)
    * @param q_i squared permutation matrix
    * @param row_start_inclusive inclusive index of start of slice
    * @param row_end_exclusive exclusive index of end of slice
    * @return mapping of col_{i+1} -> col_{i} and p_{i+1}, mapping of row_{i+1} -> row_{i} implicit via offset start_inclusive
    */
    void get_horizontal_slice(AbstractPermutation *q_i, int row_start_inclusive, int row_end_exclusive,
                              int *cur_col_to_prev_col_mapping, AbstractPermutation *succ_pi) {

        auto ptr = 0;
        for (int col = 0; col < q_i->col_size; ++col) {
            auto row = q_i->get_row_by_col(col);
            if (row >= row_start_inclusive & row < row_end_exclusive) {
                cur_col_to_prev_col_mapping[ptr] = col;
                succ_pi->set_point(row - row_start_inclusive, ptr);

                ptr++;
            }
        }
    }

    /**
     * also set -1 to no point
     * @param shrinked
     * @param row_mapper
     * @param col_mapper
     * @param flattened
     */
    inline void
    inverse_mapping(AbstractPermutation *shrinked, const int *row_mapper, const int *col_mapper,
                    AbstractPermutation *flattened) {
        //could be parallelized
        flattened->unset_all();

//#pragma omp parallel for
        for (int cur_col = 0; cur_col < shrinked->col_size; ++cur_col) {
            auto old_col = col_mapper[cur_col]; /// consecutive access

            auto cur_row = shrinked->get_row_by_col(cur_col); // consecutive access
            auto old_row = row_mapper[cur_row]; // non consecutive access
            flattened->set_point(old_row, old_col);
        }
    }


    AbstractPermutation *steady_ant(AbstractPermutation *p, AbstractPermutation *q) {
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

        auto p_lo_row_mapper = new int[spliter];
        auto p_lo = new Permutation(spliter, spliter);

        get_vertical_slice(p, 0, spliter, p_lo_row_mapper, p_lo);

        auto p_hi_row_mapper = new int[p->col_size - spliter];
        auto p_hi = new Permutation(p->col_size - spliter, p->col_size - spliter);

        get_vertical_slice(p, spliter, p->col_size, p_hi_row_mapper, p_hi);


        auto q_lo_col_mapper = new int[spliter];
        auto q_lo = new Permutation(spliter, spliter);
        get_horizontal_slice(q, 0, spliter, q_lo_col_mapper, q_lo);

        auto r_lo = p; // reuse since we no need to use p further


        auto product = steady_ant(p_lo, q_lo);


        inverse_mapping(product, p_lo_row_mapper, q_lo_col_mapper, r_lo);
//        delete product;


        auto q_hi_col_mapper = new int[p->col_size - spliter];
        auto q_hi = new Permutation(p->col_size - spliter, p->col_size - spliter);
        get_horizontal_slice(q, spliter, q->row_size, q_hi_col_mapper, q_hi);
        auto r_hi = q; // reuse since we no need to use p further
        product = steady_ant(p_hi, q_hi);


        inverse_mapping(product, p_hi_row_mapper, q_hi_col_mapper, r_hi);
//        delete product;

        //ant passage
        auto end_row = -1;
        auto end_col = r_lo->col_size + 1;
        auto cur_row = r_hi->row_size;
        auto cur_col = -1;

        auto rhi = 0;
        auto rlo = 0;

        std::vector<int> good_row_pos;
        std::vector<int> good_col_pos;
//        good_row_pos.reserve(n / 2);
//        good_col_pos.reserve(n / 2);

        bool is_went_right = false; // went up
        while (true) {

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
                std::cout << "Impossible" << std::endl;
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

        // merge r_lo to r_hi
        for (int i = 0; i < n; ++i) {
            auto col = r_lo->get_col_by_row(i);
            if (col == NOPOINT) continue;
            r_hi->set_point(i, col);
        }

        // add good points
        for (int i = 0; i < good_col_pos.size(); ++i) {
            auto col = good_col_pos[i];
            auto row = good_row_pos[i];
            r_hi->set_point(row, col);
        }

        return r_hi;
    }


    void steady_ant_with_precalc_and_memory(
            AbstractPermutation *p, AbstractPermutation *q, int *memory_block_matrices, int *free_space_matrices,
            int *memory_block_indices, PrecalcMap &map, int total_memory) {
        auto n = p->row_size;

        if (n <= map.size()) {
            auto precalced = PermutationPreAllocated(n, n, memory_block_matrices, memory_block_matrices + n,
                                                     map[n][std::hash<AbstractPermutation>()(
                                                             *p)][std::hash<AbstractPermutation>()(*q)]);
            return;
        }

        int spliter = n / 2;


        auto p_lo_row_mapper = memory_block_indices;
//        auto p_lo_row_mapper = new int[spliter];
        auto p_lo = PermutationPreAllocated(spliter, spliter, free_space_matrices, free_space_matrices + spliter);
        get_vertical_slice(p, 0, spliter, p_lo_row_mapper, &p_lo);


        auto q_lo_col_mapper = memory_block_indices + spliter;
//        auto q_lo_col_mapper = new int[spliter];
        auto q_lo = PermutationPreAllocated(spliter, spliter, free_space_matrices + 2 * spliter,
                                            free_space_matrices + 3 * spliter);
        get_horizontal_slice(q, 0, spliter, q_lo_col_mapper, &q_lo);


        auto p_hi_row_mapper = memory_block_indices + 2 * spliter;
//        auto p_hi_row_mapper = new int[n-spliter];
        auto p_hi = PermutationPreAllocated(n - spliter, n - spliter, free_space_matrices + 4 * spliter,
                                            free_space_matrices + 4 * spliter + (n - spliter));
        get_vertical_slice(p, spliter, n, p_hi_row_mapper, &p_hi);


        auto q_hi_col_mapper = memory_block_indices + 2 * spliter + (n - spliter);
//        auto q_hi_col_mapper = new int[n-spliter];
        auto q_hi = PermutationPreAllocated(n - spliter, n - spliter,
                                            free_space_matrices + 4 * spliter + 2 * (n - spliter),
                                            free_space_matrices + 4 * spliter + 3 * (n - spliter));
        get_horizontal_slice(q, spliter, q->row_size, q_hi_col_mapper, &q_hi);

        // now we have small matrices in free space, and p,q may be overwritten
//        free_space_mappings + 2*n

        // hack
        auto r_lo = p;
        auto r_hi = q;

        int on_parts = (total_memory - 2 * n) / 2;
        //
        // Permutation:
        // pairs [(row_1,col_1),....,...(row_n,col_n)]
        //


//#pragma omp parallel
//        {
//#pragma omp single
//            {
//#pragma omp task
        steady_ant_with_precalc_and_memory(&p_lo, &q_lo, free_space_matrices, memory_block_matrices,
                                           memory_block_indices + 2 * n, map, on_parts);
//#pragma omp task
        steady_ant_with_precalc_and_memory(&p_hi, &q_hi, free_space_matrices + 4 * spliter,
                                           memory_block_matrices + 4 * spliter,
                                           memory_block_indices + 2 * n + on_parts, map, on_parts);
//#pragma omp taskwait
//            }
//        };


        inverse_mapping(&p_lo, p_lo_row_mapper, q_lo_col_mapper, r_lo);
        inverse_mapping(&p_hi, p_hi_row_mapper, q_hi_col_mapper, r_hi);

        //ant passage
        auto end_row = -1;
        auto end_col = r_lo->col_size + 1;
        auto cur_row = r_hi->row_size;
        auto cur_col = -1;

        auto rhi = 0;
        auto rlo = 0;

        std::vector<int> good_row_pos;
        std::vector<int> good_col_pos;
//        good_row_pos.reserve(n / 2);
//        good_col_pos.reserve(n / 2);

        bool is_went_right = false; // went up
        while (true) {

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

        // merge r_hi to r_lo
        for (int i = 0; i < n; ++i) {
            auto col = r_hi->get_col_by_row(i);
            if (col == NOPOINT) continue;
            r_lo->set_point(i, col);
        }

        // add good points
        for (int i = 0; i < good_col_pos.size(); ++i) {
            auto col = good_col_pos[i];
            auto row = good_row_pos[i];
            r_lo->set_point(row, col);
        }

        // new matrix in r_lo

        return; //r_lo;
    }


    void
    staggered_sticky_multiplication(AbstractPermutation *p, AbstractPermutation *q, int k,
                                    PrecalcMap &map, AbstractPermutation *product) {
        if (k == p->col_size && k == q->row_size) {
            std::cout << "This function should not be called for this case, handled separately";
//            return;
        }

        if (k == 0) {
            //|q..|
            //|. p|
            for (int i = 0; i < p->row_size; ++i) {
                auto col = p->get_col_by_row(i);
                if (col != NOPOINT) product->set_point(i + q->row_size, col + q->col_size);
            }
            for (int i = 0; i < q->row_size; ++i) {
                auto col = q->get_col_by_row(i);
                if (col != NOPOINT) product->set_point(i, col);
            }
            return;
        }


        int nearest_2_degree = pow(2, int(ceil(log2(2 * k))));
        int total = int(log2(nearest_2_degree)) * nearest_2_degree;

        auto memory_block = new int[k * 8 + int(log2(nearest_2_degree)) * nearest_2_degree];
        auto mapping_row = new int[k];
        auto mapping_col = new int[k];


        auto p_red = PermutationPreAllocated(k, k, memory_block, memory_block + k);
        auto q_red = PermutationPreAllocated(k, k, memory_block + 2 * k, memory_block + 3 * k);


        // take first k columns from P and last k rows from Q, multiply and to bottom left corner of extended matrix
        get_vertical_slice(p, 0, k, mapping_row, &p_red);
        get_horizontal_slice(q, q->row_size - k, q->row_size, mapping_col, &q_red);

        steady_ant_with_precalc_and_memory(&p_red, &q_red, memory_block, memory_block + 4 * k,
                                           memory_block + 8 * k, map, total);


        // res in p_red
        for (int i = 0; i < p_red.row_size; i++) {
            auto old_col = p_red.get_col_by_row(i);
            auto cur_col = mapping_col[old_col];
            auto cur_row = mapping_row[i];
            product->set_point(q->row_size - k + cur_row, cur_col);
        }

        for (int i = 0; i < q->row_size - k; i++) {
            auto col = q->get_col_by_row(i);
            if (col != NOPOINT) product->set_point(i, col);
        }


        for (int j = k; j < p->col_size; j++) {
            auto row = p->get_row_by_col(j);
            if (row != NOPOINT) product->set_point(row + q->row_size - k, j + q->col_size - k);
        }


        delete[] memory_block;
        delete[] mapping_col;
        delete[] mapping_row;


    }

}


namespace semi_local_lcs {
    using namespace steady_ant;

    /**
    * see theorem 5.21
     * Allows get P_{a,b} when you have P_{b,a}
    */
    void fill_permutation_ba(AbstractPermutation *ab, AbstractPermutation *ba, int m, int n) {
        ba->unset_all();
        for (int i = 0; i < ab->row_size; ++i) {
            auto col = ab->get_col_by_row(i);
            if (col != NOPOINT) ba->set_point(n + m - 1 - i, m + n - 1 - col);
        }
    }

    AbstractPermutation *get_semi_local_kernel(int *a, int m, int *b, int n, PrecalcMap &map) {
        if (m == 1 && n == 1) {
            auto p = new Permutation(2, 2);
            if (a[0] == b[0]) {
                p->set_point(0, 0);
                p->set_point(1, 1);
            } else {
                p->set_point(0, 1);
                p->set_point(1, 0);
            }
            return p;

        }

        if (n > m) {
            auto n1 = n / 2;
            auto b1 = b;
            auto b2 = b + n1;

            auto subtree_l = get_semi_local_kernel(b1, n1,a, m,  map);
            auto subtree_r = get_semi_local_kernel(b2, n - n1,a, m,  map);
            auto product = new Permutation(subtree_l->row_size + subtree_r->row_size - m,
                                           subtree_l->col_size + subtree_r->col_size - m);

            auto product_t = new Permutation(subtree_l->col_size + subtree_r->col_size - m,
                                             subtree_l->row_size + subtree_r->row_size - m);


            staggered_sticky_multiplication(subtree_l, subtree_r, m, map, product);
            fill_permutation_ba(product, product_t, m, n);
            return product_t;
        } else {
            auto m1 = m / 2;
            auto a1 = a;
            auto a2 = a + m1;

            auto subtree_l = get_semi_local_kernel(a1, m1, b, n, map);
            auto subtree_r = get_semi_local_kernel(a2, m - m1, b, n, map);
            auto product = new Permutation(subtree_l->row_size + subtree_r->row_size - n,
                                           subtree_l->col_size + subtree_r->col_size - n);
            staggered_sticky_multiplication(subtree_l, subtree_r, n, map, product);

            return product;


        }

    }


};


#endif //CPU_STEADY_ANT_H
