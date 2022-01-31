#pragma once

#include <cmath>
#include <climits>
#include "Matrices.h"
#include "DominanceSum.h"

namespace common {

    /**
    * Data structure that contains pre calced product for pair of matrices i.e
    * For matrices p,q of size n it contains its product
    */
    typedef std::unordered_map<int64_t, std::unordered_map<int64_t, std::unordered_map<int64_t, int64_t>>> PrecalcMap;

    class BraidMultiplicationStrategy {
    public:
        virtual void multiply(const Permutation &m, const Permutation &n, Permutation &res) = 0;
    };

    /**
    * A naive implementation of unit Monge distance product.
    *  First, for each matrix a dominance sum matrices are constructed.
    * Then O(n^3) time complexity  distance product of explicit matrices is performed
    */
    class NaiveBraidMultiplication : public BraidMultiplicationStrategy {
    public:
        void multiply(const Permutation &p, const Permutation &q, Permutation &res) override {

            auto topRightSummator = [](int cur_row, int cur_col, int row_bound, int col_bound) {
                return cur_row >= row_bound && cur_col < col_bound;
            };

            auto dominanceM = Matrix(p.rows + 1, p.cols + 1);
            getDominanceSum(p, topRightSummator, dominanceM);
            auto dominanceN = Matrix(q.rows + 1, q.cols + 1);
            getDominanceSum(q, topRightSummator, dominanceN);

            auto rowSize = (p.rows + 1);
            auto colSize = (q.cols + 1);
            auto dominanceC = Matrix(rowSize, colSize);

            for (int i = 0; i < rowSize; ++i) {
                for (int k = 0; k < colSize; ++k) {
                    auto tmp = INT_MAX;
                    for (int j = 0; j < (p.cols + 1); ++j)
                        tmp = std::min(dominanceM.getElementAt(i, j) + dominanceN.getElementAt(j, k), tmp);
                    dominanceC.set_element_at(i, k, tmp);
                }
            }

            for (int i = 0; i < p.rows; ++i) {
                for (int j = 0; j < q.cols; ++j) {
                    auto crossDiff =
                            (dominanceC.getElementAt(i, j + 1) + dominanceC.getElementAt(i + 1, j)) -
                            (dominanceC.getElementAt(i, j) + dominanceC.getElementAt(i + 1, j + 1));
                    if (crossDiff == 1) res.set(i, j);
                }
            }
        }
    };


    /**
    * Non optimized version of steady ant algorithm that allocated new memory on each recursion step and
    * uses  precomputed if have
    * @param p
    * @param q
    * @param map
    * @return
    */
    class StickyBraidMultiplication : public BraidMultiplicationStrategy {
    public:

        StickyBraidMultiplication(PrecalcMap precalcMap) : BraidMultiplicationStrategy(), map(std::move(precalcMap)) {}


        /**
        * For each size n (1,2,3...max_size) function compute
        * for all possible permutations p of size n  (n!) and all possible permutations of size n (n!)
        * it computes its distance product.
        * For example, map[3][hash_p][hash_q] will store product of permutation matrices p and q of size 3x3
        * @param precalcMap
        * @param max_size
        */
        static void buildPrecalcMap(BraidMultiplicationStrategy *solver, PrecalcMap &precalcMap, int max_size) {
            int p_arr[max_size];
            int q_arr[max_size];
            auto empty_map = PrecalcMap();
            for (int size = 1; size < max_size + 1; size++) {
                for (int i = 0; i < size; ++i) p_arr[i] = i;

                do {
                    for (int i = 0; i < size; ++i) q_arr[i] = i;
                    do {
                        auto p = Permutation(size, size);
                        for (int i = 0; i < size; ++i) p.set(i, p_arr[i]);

                        auto q = Permutation(size, size);
                        for (int i = 0; i < size; ++i) q.set(i, q_arr[i]);

                        auto hash_p = encode(p);
                        auto hash_q = encode(q);
                        Permutation r(size, size);
                        solver->multiply(p, q, r);
                        auto points = std::vector<std::pair<int, int>>();
                        if (precalcMap[size][hash_p].count(hash_q) > 0) {
                            std::cout << " Some error";
                            return;
                        }
                        r.toPoints(points);
                        precalcMap[size][hash_p][hash_q] = encode(r);
                    } while (std::next_permutation(q_arr, q_arr + size));

                } while (std::next_permutation(p_arr, p_arr + size));
            }
        }


    protected:
        PrecalcMap map;


        /**
        * We encode permutation matrices naively. It only works till specified size (around 7-8).
        * Idea is simple, we embedded cols position in each row of non-zero entires in 32(64) bit word
        */
        static int64_t encode(const Permutation &k) {
            int64_t sum = 0;
            auto bitsPerSymbol = int(std::ceil(log2(k.rows)));

            for (int i = 0; i < k.rows; ++i) {
                sum = (sum << bitsPerSymbol) + k.getColByRow(i);
            }
            return sum;
        }

        static void decode(int64_t encoded, Permutation &k) {
            int64_t sum = 0;
            auto bitsPerSymbol = int(std::ceil(log2(k.rows)));
            auto mask = (1 << bitsPerSymbol) - 1;

            for (int row = k.rows - 1; row >= 0; row--, encoded >>= bitsPerSymbol) {
                k.set(row, encoded & mask);
            }
        }

        /**
        * see theorem 5.21
        * Allows get P_{b,a} when you have P_{a,b}
        */
        void fillPermBA(common::Permutation &ab, common::Permutation &ba, int m, int n) {
            ba.resetAll();
            for (int i = 0; i < ab.rows; ++i) {
                auto col = ab.getColByRow(i);
                if (col != NOPOINT) ba.set(n + m - 1 - i, m + n - 1 - col);
            }
        }

        /**
        * Given squared permutation matrix p_{i} get slice p[,start_inclusive:end_exclusive] and map it to new coordinates
        * to get new squared matrix of size end_exclusive-start_inclusive and mapping of row coordinates (new to old)
        * @param pI squared permutation matrix
        * @param colStartIncl inclusive index of start of slice
        * @param colEndExcl exclusive index of end of slice
        * @return mapping of row_{i+1} -> row_{i} and p_{i+1}, mapping of col_{i+1} -> col_{i} implicit via offset start_inclusive
        */
        void getVerticalSlice(const Permutation &pI, int colStartIncl, int colEndExcl, int *curRowToPrevRowMapping, Permutation &succPi) {
            auto ptr = 0;
            for (int row = 0; row < pI.rows; ++row) {
                auto col = pI.getColByRow(row);
                if (col >= colStartIncl && col < colEndExcl) {
                    curRowToPrevRowMapping[ptr] = row;
                    succPi.set(ptr, col - colStartIncl);
                    ptr++;
                }
            }
        }


        /**
         * Given squared permutation matrix q_{i} get slice p[start_inclusive:end_exclusive,:] and map it to new coordinates
        * to get new squared matrix of size end_exclusive-start_inclusive and mapping of col coordinates (new to old)
        * @param qI squared permutation matrix
        * @param rowStartIncl inclusive index of start of slice
        * @param rowEndExcl exclusive index of end of slice
        * @return mapping of col_{i+1} -> col_{i} and p_{i+1}, mapping of row_{i+1} -> row_{i} implicit via offset start_inclusive
        */
        void getHorizontalSlice(const Permutation &qI, int rowStartIncl, int rowEndExcl, int *curColToPrevColMapping, Permutation &succPi) {
            auto ptr = 0;
            for (int col = 0; col < qI.cols; ++col) {
                auto row = qI.getRowByCol(col);
                if (row >= rowStartIncl && row < rowEndExcl) {
                    curColToPrevColMapping[ptr] = col;
                    succPi.set(row - rowStartIncl, ptr);
                    ptr++;
                }
            }
        }

        /**
         * Maps non-zero entries of shrinked n/2Xn/2 permutation matrix to nXn matrix. Aka maps positions of non-zero
         * elements of matrix shrinked to flattened matrix  according to row_mapper and col_mapper
         * @param shrinked
         * @param rowMapper
         * @param colMapper
         * @param flattened
         */
        inline void inverseMapping(const Permutation &shrinked, const int *rowMapper, const int *colMapper, Permutation &flattened) {
            //could be parallelized
            flattened.resetAll();

            //#pragma omp parallel for
            for (int curCol = 0; curCol < shrinked.cols; ++curCol) {
                auto oldCol = colMapper[curCol]; /// consecutive access
                auto curRow = shrinked.getRowByCol(curCol); // consecutive access
                auto oldRow = rowMapper[curRow]; // non consecutive access
                flattened.set(oldRow, oldCol);
            }
        }

        /**
        * Maps non-zero entries of shrinked n/2Xn/2 permutation matrix to nXn matrix. Aka maps positions of non-zero
        * elements of matrix shrinked to flattened matrix  according to row_mapper and col_mapper
        * @param shrinked
        * @param row_mapper
        * @param col_mapper
        * @param flattened
        */
        inline void
        inverseMappingWithoutPreClearing(Permutation &shrinked, const int *row_mapper,
                                         const int *col_mapper,
                                         Permutation &flattened) {

            //#pragma omp parallel for
            for (int curCol = 0; curCol < shrinked.cols; ++curCol) {
                auto oldCol = col_mapper[curCol]; /// consecutive access

                auto curRow = shrinked.getRowByCol(curCol); // consecutive access
                auto oldRow = row_mapper[curRow]; // non consecutive access
                flattened.set(oldRow, oldCol);
            }
        }


        inline void antPassage(Permutation &rLo, Permutation &rHi, int n, std::vector<int> &goodRowPos, std::vector<int> &goodColPos) {

            //ant passage
            auto endRow = -1;
            auto endCol = rLo.cols + 1;
            auto curRow = rHi.rows;
            auto curCol = -1;

            auto rhi = 0;
            auto rlo = 0;

            bool isWentRight = false; // went up
            while (true) {

                if (endCol == curCol && endRow == curRow) break;

                if (curRow == 0) break;
                //TODO is new
                if (curCol == n) break;
                //
                auto dominanceRow = curRow - 1;
                auto dominanceCol = curCol + 1;

                //prev step
                if (isWentRight) {
                    rhi = IncrementalDominanceSum::move<BOTTMO_RIGHT, RIGHT>(dominanceRow, dominanceCol - 1, rhi, rHi);
                    rlo = IncrementalDominanceSum::move<TOP_LEFT, RIGHT>(dominanceRow, dominanceCol - 1, rlo, rLo);
                } else {
                    rhi = IncrementalDominanceSum::move<BOTTMO_RIGHT, UP>(dominanceRow + 1, dominanceCol, rhi, rHi);
                    rlo = IncrementalDominanceSum::move<TOP_LEFT, UP>(dominanceRow + 1, dominanceCol, rlo, rLo);
                }

                if (rhi - rlo < 0) {
                    isWentRight = true;
                    curCol++;
                } else if (rhi - rlo == 0) {
                    isWentRight = false;
                    curRow--;
                } else {
                    std::cout << "Impissble" << std::endl;
                }

                if (dominanceCol > 0) {
                    auto delta_above_left = IncrementalDominanceSum::move<BOTTMO_RIGHT, LEFT>(dominanceRow, dominanceCol, rhi, rHi) -
                                            IncrementalDominanceSum::move<TOP_LEFT, LEFT>(dominanceRow, dominanceCol, rlo, rLo);

                    auto delta_below_right = IncrementalDominanceSum::move<BOTTMO_RIGHT, DOWN>(dominanceRow, dominanceCol, rhi, rHi) -
                                             IncrementalDominanceSum::move<TOP_LEFT, DOWN>(dominanceRow, dominanceCol, rlo, rLo);

                    if (delta_above_left < 0 && delta_below_right > 0) {
                        goodRowPos.push_back(dominanceRow);
                        goodColPos.push_back(dominanceCol - 1);
                    }
                }
            }
            // end ant passage
        }
    };

    class SimpleStickyBraidMultiplication : public StickyBraidMultiplication {

    public:
        SimpleStickyBraidMultiplication(PrecalcMap precalcMap) : StickyBraidMultiplication(std::move(precalcMap)) {}

        void multiply(const Permutation &p, const Permutation &q, Permutation &res) override {
            auto pCopy = p;
            auto qCopy = q;
            res = mul(pCopy, qCopy);
        }

    private:

        Permutation mul(Permutation &p, Permutation &q) {
            auto n = p.cols;

            if (n <= this->map.size()) {
                const auto &pts = this->map[n][encode(p)][encode(q)];
                Permutation res(n, n);
                decode(pts, res);
                return res;
            }


            if (n == 1) {
                //base case when common dimension is one
                auto row = p.getRowByCol(0);
                auto col = q.getColByRow(0);
                auto res = Permutation(p.rows, q.cols);
                if (row != NOPOINT && col != NOPOINT) res.set(row, col);
                return res;
            }

            int spliter = p.cols / 2;
            auto sz = p.cols - spliter;

            auto p_lo_row_mapper = new int[spliter];
            auto p_lo = Permutation(spliter, spliter);

            this->getVerticalSlice(p, 0, spliter, p_lo_row_mapper, p_lo);

            auto p_hi_row_mapper = new int[sz];
            auto p_hi = Permutation(sz, sz);

            this->getVerticalSlice(p, spliter, p.cols, p_hi_row_mapper, p_hi);


            auto q_lo_col_mapper = new int[spliter];
            auto q_lo = Permutation(spliter, spliter);
            this->getHorizontalSlice(q, 0, spliter, q_lo_col_mapper, q_lo);

            auto &r_lo = p; // reuse since we no need to use p further
            auto product = this->mul(p_lo, q_lo);
            this->inverseMapping(product, p_lo_row_mapper, q_lo_col_mapper, r_lo);

            auto q_hi_col_mapper = new int[sz];
            auto q_hi = Permutation(sz, sz);
            this->getHorizontalSlice(q, spliter, q.rows, q_hi_col_mapper, q_hi);
            auto r_hi = q;
            product = this->mul(p_hi, q_hi);

            this->inverseMapping(product, p_hi_row_mapper, q_hi_col_mapper, r_hi);

            delete[] p_lo_row_mapper;
            delete[] q_lo_col_mapper;
            delete[] p_hi_row_mapper;
            delete[] q_hi_col_mapper;

            std::vector<int> good_row_pos;
            std::vector<int> good_col_pos;
            this->antPassage(r_lo, r_hi, n, good_row_pos, good_col_pos);

            // merge r_lo to r_hi
            for (int i = 0; i < n; ++i) {
                auto col = r_lo.getColByRow(i);
                if (col == NOPOINT) continue;
                r_hi.set(i, col);
            }

            // add good points
            for (int i = 0; i < good_col_pos.size(); ++i) {
                auto col = good_col_pos[i];
                auto row = good_row_pos[i];
                r_hi.set(row, col);
            }

            return r_hi;
        }
    };

}