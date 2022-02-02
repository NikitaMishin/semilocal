#pragma once

#include <cmath>
#include <climits>
#include "Matrices.h"
#include "DominanceSum.h"

namespace common {

    class StaggeredStickyBraidMultiplication;
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
                    dominanceC.setElementAt(i, k, tmp);
                }
            }
            dominanceC.getCrossDiffPermutation(res);
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
        template<class T>
        static int64_t encode(const T &k) {
            int64_t sum = 0;
            auto bitsPerSymbol = int(std::ceil(log2(k.rows)));

            for (int i = 0; i < k.rows; ++i) {
                sum = (sum << bitsPerSymbol) + k.getColByRow(i);
            }
            return sum;
        }

        template<class T>
        static void decode(int64_t encoded, T &k) {
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
        template<class T>
        void fillPermBA(T &ab, T &ba, int m, int n) {
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
        template<class T>
        void getVerticalSlice(const T &pI, int colStartIncl, int colEndExcl, int *curRowToPrevRowMapping, T &succPi) {
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
        template<class T>
        void getHorizontalSlice(const T &qI, int rowStartIncl, int rowEndExcl, int *curColToPrevColMapping, T &succPi) {
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
        template<class T>
        inline void inverseMapping(const T &shrinked, const int *rowMapper, const int *colMapper, T &flattened) {
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
        template<class T>
        inline void inverseMappingWithoutPreClearing(T &shrinked, const int *row_mapper, const int *col_mapper, T &flattened) {

            //#pragma omp parallel for
            for (int curCol = 0; curCol < shrinked.cols; ++curCol) {
                auto oldCol = col_mapper[curCol]; /// consecutive access

                auto curRow = shrinked.getRowByCol(curCol); // consecutive access
                auto oldRow = row_mapper[curRow]; // non consecutive access
                flattened.set(oldRow, oldCol);
            }
        }


        template<class T>
        inline void antPassage(const T &rLo, const T &rHi, int n, std::vector<int> &goodRowPos, std::vector<int> &goodColPos) {

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


        template<class T>
        inline void restoreMatrixInRlo(T &rLo, const T &rHi, int n) {
            std::vector<int> goodRowPos;
            std::vector<int> goodColPos;
            this->antPassage(rLo, rHi, n, goodRowPos, goodColPos);

            // merge r_hi to r_lo
            for (int i = 0; i < n; ++i) {
                auto col = rHi.getColByRow(i);
                if (col == NOPOINT) continue;
                rLo.set(i, col);
            }

            // add good points
            for (int i = 0; i < goodColPos.size(); ++i) {
                auto col = goodColPos[i];
                auto row = goodRowPos[i];
                rLo.set(row, col);
            }
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

            auto pLoRowMapper = new int[spliter];
            auto pLo = Permutation(spliter, spliter);

            this->getVerticalSlice(p, 0, spliter, pLoRowMapper, pLo);

            auto pHiRowMapper = new int[sz];
            auto pHi = Permutation(sz, sz);

            this->getVerticalSlice(p, spliter, p.cols, pHiRowMapper, pHi);


            auto qLoColMapper = new int[spliter];
            auto qLo = Permutation(spliter, spliter);
            this->getHorizontalSlice(q, 0, spliter, qLoColMapper, qLo);

            auto &rLo = p; // reuse since we no need to use p further
            auto product = this->mul(pLo, qLo);
            this->inverseMapping(product, pLoRowMapper, qLoColMapper, rLo);

            auto qHiColMapper = new int[sz];
            auto qHi = Permutation(sz, sz);
            this->getHorizontalSlice(q, spliter, q.rows, qHiColMapper, qHi);
            auto rHi = q;
            product = this->mul(pHi, qHi);

            this->inverseMapping(product, pHiRowMapper, qHiColMapper, rHi);

            delete[] pLoRowMapper;
            delete[] qLoColMapper;
            delete[] pHiRowMapper;
            delete[] qHiColMapper;

            restoreMatrixInRlo(rLo, rHi, n);

            return rLo;
        }
    };


    class MemoryOptimizedStickBraidMultiplication : public StickyBraidMultiplication {
    public:
        friend StaggeredStickyBraidMultiplication;

        MemoryOptimizedStickBraidMultiplication(PrecalcMap precalcMap) : StickyBraidMultiplication(std::move(precalcMap)) {}

    protected:
        /**
        * For directly managing raw block of memories
        */
        class PermutationPreAllocated {

        public:

            int *rowToCol;
            int *colToRow;
            int rows;
            int cols;

            inline void set(int row, int col) {
                rowToCol[row] = col;
                colToRow[col] = row;
            }

            inline void reset(int row, int col) {
                if (rowToCol[col] == col) {
                    rowToCol[col] = NOPOINT;
                    colToRow[row] = NOPOINT;
                }
            }

            inline void resetAll() {
                for (int i = 0; i < rows; ++i) rowToCol[i] = NOPOINT;
                for (int i = 0; i < cols; ++i) colToRow[i] = NOPOINT;
            }


            inline int getRowByCol(int col) const { return colToRow[col]; }

            inline int getColByRow(int row) const { return rowToCol[row]; }


            void toPoints(std::vector<std::pair<int, int>> &result) const {
                for (int i = 0; i < cols; ++i) {
                    auto col = getColByRow(i);
                    if (col != NOPOINT) result.emplace_back(i, col);
                }
            }

            PermutationPreAllocated() = default;

            PermutationPreAllocated(const PermutationPreAllocated &other) {
                rows = other.rows;
                cols = other.cols;
                rowToCol = other.rowToCol;
                colToRow = other.colToRow;
            }

            PermutationPreAllocated(int row, int col, int *row_to_col, int *col_to_row) : rowToCol(row_to_col), colToRow(col_to_row), rows(row), cols(col) {};

            PermutationPreAllocated(int row, int col, int *row_to_col, int *col_to_row,
                                    std::vector<std::pair<int, int>> &points) : PermutationPreAllocated(row, col, row_to_col, col_to_row) {
                for (int i = 0; i < rows; ++i) row_to_col[i] = NOPOINT;
                for (int i = 0; i < cols; ++i) col_to_row[i] = NOPOINT;
                for (auto &point: points) {
                    row_to_col[point.first] = point.second;
                    col_to_row[point.second] = point.first;
                }
            }
        };

        template<class F, class T>
        inline void copy(const F &from, T &to) {
            for (int i = 0; i < from.rows; ++i) {
                auto col = from.getColByRow(i);
                to.set(i, col);
            }
        }

        using P = PermutationPreAllocated;
        struct MemoryBlocks {
            int *matrices;
            int *freeSpaceMatrices;
            int *indices;
        };

        bool initialPart(int n, const P &p, const P &q, MemoryBlocks mem, int *&pLoRowMapper, P &pLo, int *&qLoColMapper, P &qLo, int *&pHiRowMapper, P &pHi,
                         int *&qHiColMapper, P &qHi) {
            if (n <= map.size()) {
                auto precalced = PermutationPreAllocated(n, n, mem.matrices, mem.matrices + n);
                decode(map[n][encode(p)][encode(q)], precalced);//rewrite memory
                return true;
            }

            int spliter = n / 2;


            pLoRowMapper = mem.indices;
            pLo = PermutationPreAllocated(spliter, spliter, mem.freeSpaceMatrices, mem.freeSpaceMatrices + spliter);
            getVerticalSlice(p, 0, spliter, pLoRowMapper, pLo);

            qLoColMapper = mem.indices + spliter;
            qLo = PermutationPreAllocated(spliter, spliter, mem.freeSpaceMatrices + 2 * spliter, mem.freeSpaceMatrices + 3 * spliter);
            getHorizontalSlice(q, 0, spliter, qLoColMapper, qLo);

            pHiRowMapper = mem.indices + 2 * spliter;
            pHi = PermutationPreAllocated(n - spliter, n - spliter, mem.freeSpaceMatrices + 4 * spliter, mem.freeSpaceMatrices + 4 * spliter + (n - spliter));
            getVerticalSlice(p, spliter, n, pHiRowMapper, pHi);


            qHiColMapper = mem.indices + 2 * spliter + (n - spliter);
            qHi = PermutationPreAllocated(n - spliter, n - spliter,
                                          mem.freeSpaceMatrices + 4 * spliter + 2 * (n - spliter),
                                          mem.freeSpaceMatrices + 4 * spliter + 3 * (n - spliter));
            getHorizontalSlice(q, spliter, q.rows, qHiColMapper, qHi);
            return false;
        }


        std::pair<MemoryBlocks, int> allocate(int k, bool parallel) {
            int nearest2Degree = pow(2, int(ceil(log2(2 * k))));

            // then we need to use O(nlogn) memory if parallel for indices
            int memOnIndices = (parallel) ? int(log2(nearest2Degree)) * nearest2Degree : k * 8;
            int *memoryBlock;
            memoryBlock = new int[k * 8 + memOnIndices];
            auto freeBlock1 = memoryBlock;
            auto freeBlock2 = memoryBlock + 4 * k;
            auto freeIndicesBlock = memoryBlock + 8 * k;

            return {{.matrices=freeBlock1, .freeSpaceMatrices=freeBlock2, .indices=freeIndicesBlock}, memOnIndices};

        }
    };

    class SequentialMemoryOptimizedStickBraidMultiplication : public MemoryOptimizedStickBraidMultiplication {
    public:

        SequentialMemoryOptimizedStickBraidMultiplication(PrecalcMap precalcMap) : MemoryOptimizedStickBraidMultiplication(std::move(precalcMap)) {}

        void multiply(const Permutation &p, const Permutation &q, Permutation &res) override {
            int k = p.rows;
            auto[memory, _] = allocate(p.rows, false);

            auto pRed = PermutationPreAllocated(k, k, memory.matrices, memory.matrices + k);
            auto qRed = PermutationPreAllocated(k, k, memory.matrices + 2 * k, memory.matrices + 3 * k);
            copy(p, pRed);
            copy(q, qRed);
            steadyAntWithPrecalcAndMemory(pRed, qRed, memory);
            res = Permutation(p.rows, q.cols);
            copy(pRed, res);
            delete memory.matrices;
        }

    private:

        void
        steadyAntWithPrecalcAndMemory(PermutationPreAllocated &p, PermutationPreAllocated &q, MemoryBlocks mem) {
            auto n = p.rows;
            int *pLoRowMapper, *qLoColMapper, *pHiRowMapper, *qHiColMapper;
            P pLo, qLo, pHi, qHi;
            auto baseCase = initialPart(n, p, q, mem, pLoRowMapper, pLo, qLoColMapper, qLo, pHiRowMapper, pHi, qHiColMapper, qHi);
            if (baseCase) return;
            int spliter = n / 2;

            // hack
            auto rLo = p;
            auto rHi = q;

            MemoryBlocks loMem{.matrices = mem.freeSpaceMatrices, .freeSpaceMatrices = mem.matrices, .indices = mem.indices + 2 * n};
            steadyAntWithPrecalcAndMemory(pLo, qLo, loMem);
            MemoryBlocks hiMem{.matrices = mem.freeSpaceMatrices + 4 * spliter,
                    .freeSpaceMatrices = mem.matrices + 4 * spliter,
                    .indices = mem.indices + 2 * n};
            steadyAntWithPrecalcAndMemory(pHi, qHi, hiMem);

            inverseMapping(pLo, pLoRowMapper, qLoColMapper, rLo);
            inverseMapping(pHi, pHiRowMapper, qHiColMapper, rHi);

            restoreMatrixInRlo(rLo, rHi, n);
        }
    };


    class OpenMPStickyBraid : public MemoryOptimizedStickBraidMultiplication {

    public:
        OpenMPStickyBraid(int nestedParallRegions, PrecalcMap precalcMap) : MemoryOptimizedStickBraidMultiplication(std::move(precalcMap)),
                                                                            parallelizationFactor(nestedParallRegions) {}

        void multiply(const Permutation &p, const Permutation &q, Permutation &res) override {
            int k = p.rows;
            auto[memory, total] = allocate(p.rows, true);

            auto p_red = PermutationPreAllocated(k, k, memory.matrices, memory.matrices + k);
            auto q_red = PermutationPreAllocated(k, k, memory.matrices + 2 * k, memory.matrices + 3 * k);
            copy(p, p_red);
            copy(q, q_red);
            parall(p_red, q_red, memory, total, parallelizationFactor);
            res = Permutation(p.rows, q.cols);
            copy(p_red, res);
            delete memory.matrices;
        }

    private:
        int parallelizationFactor;

        void parall(P &p, P &q, MemoryBlocks mem, int totalMemory, int nestedParallRegions = 2) {
            auto n = p.rows;
            int *pLoRowMapper, *qLoColMapper, *pHiRowMapper, *qHiColMapper;
            P pLo, qLo, pHi, qHi;
            auto baseCase = initialPart(n, p, q, mem, pLoRowMapper, pLo, qLoColMapper, qLo, pHiRowMapper, pHi, qHiColMapper, qHi);
            if (baseCase) return;
            int spliter = n / 2;

            // hack
            auto rLo = p;
            auto rHi = q;

            int on_parts = (totalMemory - 2 * n) / 2;
            if (nestedParallRegions > 0) {

#pragma omp parallel num_threads(2)
                {
#pragma omp single nowait
                    {

#pragma omp task depend(in: on_parts)
                        parall(pLo, qLo,
                               MemoryBlocks{.matrices = mem.freeSpaceMatrices, .freeSpaceMatrices = mem.matrices, .indices=mem.indices + 2 * n}, on_parts,
                               nestedParallRegions - 1);
#pragma omp task depend(in: on_parts)
                        parall(pHi, qHi, MemoryBlocks{.matrices = mem.freeSpaceMatrices + 4 * spliter,
                                       .freeSpaceMatrices = mem.matrices + 4 * spliter, .indices=mem.indices + 2 * n + on_parts}, on_parts,
                               nestedParallRegions - 1);

#pragma omp task depend(out: on_parts)
                        {
                            inverseMapping(pLo, pLoRowMapper, qLoColMapper, rLo);
                            inverseMapping(pHi, pHiRowMapper, qHiColMapper, rHi);
                            restoreMatrixInRlo(rLo, rHi, n);
                        }
//#pragma omp taskwait
                    };
                }

            } else {
                parall(pLo, qLo,
                       MemoryBlocks{.matrices = mem.freeSpaceMatrices, .freeSpaceMatrices = mem.matrices, .indices=mem.indices + 2 * n}, on_parts,
                       nestedParallRegions);
                parall(pHi, qHi, MemoryBlocks{.matrices = mem.freeSpaceMatrices + 4 * spliter, .freeSpaceMatrices = mem.matrices + 4 * spliter,
                        .indices=mem.indices + 2 * n + on_parts}, on_parts, nestedParallRegions);
                inverseMapping(pLo, pLoRowMapper, qLoColMapper, rLo);
                inverseMapping(pHi, pHiRowMapper, qHiColMapper, rHi);
                restoreMatrixInRlo(rLo, rHi, n);

            }
        }


    };

    class StaggeredStickyBraidMultiplication {
    public:
        StaggeredStickyBraidMultiplication(MemoryOptimizedStickBraidMultiplication *stickyBraidSolver) : solver(stickyBraidSolver) {}

        template<bool RowGlue>
        void staggeredStickyMultiplication(const Permutation &p, Permutation &q, int k, Permutation &product) {
            if (k == p.rows && k == q.cols) {
                std::cout << "This function should not be called for this case, handled separately";
                return;
            }
            //TODO dont forget to init
            //TODO create callback with common fillment

            if (k == 0) {
                for (int i = 0; i < p.rows; ++i) {
                    auto col = p.getColByRow(i);
                    if constexpr(RowGlue) {
                        if (col != NOPOINT) product.set(i + q.rows, col + q.cols);
                    } else {
                        if (col != NOPOINT) product.set(i, col);
                    }
                }
                for (int i = 0; i < q.rows; ++i) {
                    auto col = q.getColByRow(i);
                    if constexpr(RowGlue) {
                        if (col != NOPOINT) product.set(i, col);
                    } else {
                        if (col != NOPOINT) product.set(i + p.rows, col + p.cols);
                    }
                }
                return;
            }

            auto pRed = Permutation(k, k);
            auto qRed = Permutation(k, k);
            auto mapping_row = new int[k];
            auto mapping_col = new int[k];

            if constexpr(RowGlue) {
                // take first k columns from P and last k rows from Q, multiply and to bottom left corner of extended matrix
                solver->getVerticalSlice(p, 0, k, mapping_row, pRed);
                solver->getHorizontalSlice(q, q.rows - k, q.rows, mapping_col, qRed);
            } else {
                // take last k columns from P and first k rows from Q, multiply
                solver->getVerticalSlice(p, p.rows - k, p.rows, mapping_row, pRed);
                solver->getHorizontalSlice(q, 0, k, mapping_col, qRed);
            }
            Permutation res;
            solver->multiply(pRed, qRed, res);

            for (int i = 0; i < res.rows; i++) {
                auto old_col = res.getColByRow(i);
                auto cur_col = mapping_col[old_col];
                auto cur_row = mapping_row[i];
                if (RowGlue) {
                    product.set(q.rows - k + cur_row, cur_col);
                } else {
                    product.set(cur_row, cur_col + p.cols - k);
                }
            }

            if constexpr (RowGlue) {
                for (int j = k; j < p.cols; j++) {
                    auto row = p.getRowByCol(j);
                    if (row != NOPOINT) product.set(row + q.rows - k, j + q.cols - k);
                }
                for (int i = 0; i < q.rows - k; i++) {
                    auto col = q.getColByRow(i);
                    if (col != NOPOINT) product.set(i, col);
                }
            } else {
                for (int j = 0; j < p.cols - k; j++) {
                    auto row = p.getRowByCol(j);
                    if (row != NOPOINT) product.set(row, j);
                }
                for (int i = k; i < q.rows; i++) {
                    auto col = q.getColByRow(i);
                    if (col != NOPOINT) product.set(i - k + p.rows, p.cols - k + col);//here
                }
            }
            delete[] mapping_col;
            delete[] mapping_row;
        }

        void glueingPartToWhole(Permutation &whole, Permutation &part, int offsetL, int offset_r, Permutation &product) {
            if (part.rows != (whole.rows - offsetL - offset_r)) {
                throw std::runtime_error("Dimensions not match");
            }
            auto k = whole.rows - offset_r - offsetL;

            // first offset_l strands goes inact
            for (int col = 0; col < offsetL; ++col) {
                auto row = whole.getColByRow(col);
                product.set(row, col);
            }
            // last offset_r strands goes inact
            for (int col = whole.rows - offset_r; col < whole.rows; ++col) {
                auto row = whole.getRowByCol(col);
                product.set(row, col);
            }

            auto wholeRed = Permutation(k, k);
            auto partCopy = Permutation(k, k);

            auto mappingRow = new int[k];
            // get strands that insersects  with part braid  (its strands that hit border with part braid)
            solver->getVerticalSlice(whole, offsetL, k + offsetL, mappingRow, wholeRed);
            solver->copy(part, partCopy);

            Permutation res;
            solver->multiply(wholeRed, partCopy, res);
            for (int i = 0; i < res.rows; i++) {
                auto oldCol = res.getColByRow(i);
                auto curCol = oldCol;
                auto curRow = mappingRow[i];
                product.set(curRow, curCol + offsetL);
            }

            delete[] mappingRow;
        };

    private:
        MemoryOptimizedStickBraidMultiplication *solver;
    };


}