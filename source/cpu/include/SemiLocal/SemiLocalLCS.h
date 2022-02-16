#pragma once

#include "../../common/include/Permutations/MongeMultuplicationStrategy.h"

namespace semi_local {

    template<class Input>
    class SemiLocalStrategy {
    public:
        virtual void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept;
    };

    template<class Input>
    class DummyPrefixLCSStrategy : public SemiLocalStrategy<Input> {
    public:
        DummyPrefixLCSStrategy(Input wildCardSymbol) : matchAny(wildCardSymbol) {}

        void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept override {
            auto customLCS = [a, aSize, this](Input *bSymbol, int ii, int jj) {
                auto m = aSize + 1;
                auto bExt = bSymbol + ii;
                auto n = jj - ii + 1;
                int arr[m][n];
                for (auto i = 0; i < m; i++) arr[i][0] = 0;
                for (auto j = 0; j < n; j++) arr[0][j] = 0;

                for (int i = 1; i < m; ++i) {
                    for (int j = 1; j < n; ++j) {
                        arr[i][j] = std::max(
                                std::max(arr[i - 1][j], arr[i][j - 1]),
                                (a[i - 1] == bExt[j - 1] || (bExt[j - 1] == matchAny)) ? arr[i - 1][j - 1] + 1 : arr[i - 1][j - 1]);
                    }
                }
                return arr[m - 1][n - 1];
            };

            auto bExtended = new Input[aSize + bSize + aSize];
            for (int i = 0; i < bSize; ++i) bExtended[aSize + i] = b[i];
            for (int i = 0; i < aSize; ++i) {
                bExtended[i] = matchAny;
                bExtended[bSize + aSize + i] = matchAny;
            }

            auto matrix = common::Matrix(aSize + bSize + 1, aSize + bSize + 1);
            int m = aSize;
            int n = bSize;
            for (int i = 0; i < matrix.rows; ++i) {
                for (int j = 0; j < matrix.cols; j++) {
                    if (j <= i - m) {
                        matrix.setElementAt(i, j, -(j - i + m));
                    } else {
                        matrix.setElementAt(i, j, -customLCS(bExtended, i, j + m));
                    }
                }
            }
            delete[] bExtended;
            kernel = common::Permutation(m + n, m + n);
            matrix.getCrossDiffPermutation(kernel);
        }

    private:
        Input matchAny;
    };

    template<class Input, bool WithIf>
    class SimpleIterativeCombing : public SemiLocalStrategy<Input> {
    public:
        void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept override {
            auto m = aSize;
            auto n = bSize;

            auto top_strands = new Input[n];
            auto left_strands = new Input[m];

            // init phase
            for (int i = 0; i < m; ++i) left_strands[i] = i;

            for (int i = 0; i < n; ++i) top_strands[i] = i + m;


            for (int i = 0; i < m; ++i) {
                auto left_edge = m - 1 - i;
                auto left_strand = left_strands[left_edge];
                auto a_symbol = a[i];
                int right_strand;
                for (int j = 0; j < n - 1; ++j) {
                    right_strand = top_strands[j];
                    auto r = a_symbol == b[j] || (left_strand > right_strand);

                    if constexpr(WithIf) {
                        if (r) {
                            top_strands[j] = left_strand;
                            left_strand = right_strand;
                        }
                    } else {

                        top_strands[j] = (right_strand & (r - 1)) | ((-r) & left_strand);
                        left_strand = (left_strand & (r - 1)) | ((-r) & right_strand);

                    }

                }

                right_strand = top_strands[n - 1];
                auto r = a_symbol == b[n - 1] || (left_strand > right_strand);
                left_strands[left_edge] = (left_strand & (r - 1)) | ((-r) & right_strand);

                if constexpr(WithIf) {
                    if (r) top_strands[n - 1] = left_strand;
                } else {
                    top_strands[n - 1] = (right_strand & (r - 1)) | ((-r) & left_strand);
                }

            }


            kernel = common::Permutation(n + m, n + m);
            // permutation construction phase
            for (int l = 0; l < m; l++) kernel.set(left_strands[l], n + l);

            for (int r = m; r < m + n; r++) kernel.set(top_strands[r - m], r - m);

            delete[] left_strands;
            delete[] top_strands;
        }
    };

    template<class Input, bool WithIf>
    class AntidiagonalIterativeCombing : public SemiLocalStrategy<Input> {
    public:


    protected:

        /**
         *
         * @tparam Input
         * @tparam WithIf weather or not to use approach with if rather then  the branchless one @param strandMap
         * @param a
         * @param b
         * @param upperBound
         * @param leftEdge
         * @param topEdge
         * @param offsetA
         * @param offsetB
         */
        template<bool Parallel>
        inline void
        antiDiagonalComputation(Input *strandMap, const Input *a, const Input *b, int upperBound, int leftEdge,
                                int topEdge, int offsetA, int offsetB) noexcept {

#pragma omp for simd schedule(static) aligned(a, b, strandMap:sizeof(Input)*8) nowait
            for (int k = 0; k < upperBound; ++k) {

                auto left_strand = strandMap[leftEdge + k];
                auto right_strand = strandMap[topEdge + k];

                auto r = (a[offsetA + k] == b[offsetB + k]) || (left_strand > right_strand);

                if constexpr(WithIf) {
                    if (r) {
                        strandMap[topEdge + k] = left_strand;
                        strandMap[leftEdge + k] = right_strand;
                    }
                } else {
                    auto r_minus = (r - 1);
                    auto minus_r = -r;
                    auto l_new = (left_strand & r_minus) | (minus_r & right_strand);
                    auto r_new = (right_strand & r_minus) | (minus_r & left_strand);

                    strandMap[leftEdge + k] = l_new;
                    strandMap[topEdge + k] = r_new;
                }
            }

            if constexpr(Parallel) {
#pragma omp barrier
            }
        }


        inline void initialization(Input *strandMap, int m, int n) noexcept {
#pragma omp for simd schedule(static)
            for (int k = 0; k < m; ++k) {
                strandMap[k] = k;
            }

#pragma omp for simd schedule(static)
            for (int l = 0; l < n; ++l) {
                strandMap[l + m] = l + m;
            }

        }

        inline void
        constructPermutation(common::Permutation &matrix, Input *strand_map, bool is_reverse, int m, int n) noexcept {
            if (!is_reverse) {
#pragma omp for simd schedule(static)
                for (int r = m; r < m + n; r++) {
                    matrix.set(strand_map[r], r - m);
                }
#pragma omp for simd schedule(static)
                for (int l = 0; l < m; l++) {
                    matrix.set(strand_map[l], n + l);
                }


            } else {

#pragma omp for simd schedule(static)
                for (int r = m; r < m + n; r++) {
                    matrix.set(n + m - 1 - strand_map[r], n + m - 1 - (r - m));
                }
#pragma omp for simd schedule(static)
                for (int l = 0; l < m; l++) {
                    matrix.set(n + m - 1 - strand_map[l], n + m - 1 - (n + l));
                }

            }

        }

        inline void fillAReverse(const Input *a, Input *aReverse, int m) noexcept {
#pragma omp  for simd schedule(static)
            for (int i = 0; i < m; ++i) {
                aReverse[i] = a[m - 1 - i];
            }
        }
    };

    template<class Input, bool WithIf>
    class OpenMPIterativeCombing : public AntidiagonalIterativeCombing<Input, WithIf> {
    public:

        OpenMPIterativeCombing(int numThreads) : AntidiagonalIterativeCombing<Input, WithIf>(), threadsNum(numThreads) {}


        void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept override {
            if (threadsNum > 1) {
                stickyBraidMpi<true>(kernel, a, aSize, b, bSize);
            } else {
                stickyBraidMpi<false>(kernel, a, aSize, b, bSize);
            }
        }

    private:
        int threadsNum;

        template<bool Parallel>
        void stickyBraidMpi(common::Permutation &matrix, const Input *a, int aSize, const Input *b, int bSize, bool isReverse = false) noexcept {

            if (aSize > bSize) {
                stickyBraidMpi<Parallel>(matrix, b, bSize, a, aSize, !isReverse);
                return;
            }

            auto m = aSize;
            auto n = bSize;
            matrix = common::Permutation(m + n, m + n);

            auto size = m + n;
            auto *strand_map = new Input[size];

            auto numDiag = m + n - 1;
            auto totalSameLengthDiag = numDiag - (m - 1) - (m - 1);
            auto *aReverse = new Input[m];


#pragma omp parallel num_threads(threadsNum)  default(none) shared(aReverse, a, b, isReverse, strand_map, matrix, totalSameLengthDiag, size, m, n)
            {
                int leftEdge, topEdge;
                //    init phase
                this->initialization(strand_map, m, n);
                this->fillAReverse(a, aReverse, m);

                //    phase one
                topEdge = m;
                leftEdge = m - 1;
                for (int curDiagLen = 0; curDiagLen < m - 1; ++curDiagLen) {
                    this->template antiDiagonalComputation<Parallel>(strand_map, aReverse, b, curDiagLen + 1, leftEdge, topEdge, leftEdge, 0);
                    leftEdge--;
                }

                //phase 2
                topEdge = m;
                for (int j = 0; j < totalSameLengthDiag; ++j) {
                    this->template antiDiagonalComputation<Parallel>(strand_map, aReverse, b, m, 0, topEdge, 0, j);
                    topEdge++;
                }

                //// phase 3
                auto startJ = totalSameLengthDiag;
                topEdge = startJ + m;
                for (int diagLen = m - 2; diagLen >= 0; --diagLen, startJ++) {
                    this->template antiDiagonalComputation<Parallel>(strand_map, aReverse, b, diagLen + 1, 0, topEdge, 0, startJ);
                    topEdge++;
                }

                this->constructPermutation(matrix, strand_map, isReverse, m, n);
            }

            delete[] aReverse;
            delete[] strand_map;
        }

    };

    template<class Input, bool WithIf>
    class OpenMPIterativeCombingCombined : public AntidiagonalIterativeCombing<Input, WithIf> {

    public:
        OpenMPIterativeCombingCombined(int numThreads, common::StaggeredStickyBraidMultiplication *staggeredStickyBraidMultiplicator) :
                AntidiagonalIterativeCombing<Input, WithIf>(), threadsNum(numThreads), staggeredStickyBraidMultiplication(staggeredStickyBraidMultiplicator) {}

        void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept override {
            if (threadsNum > 1) {
                firstAndThirdPhaseCombined<true>(kernel, a, aSize, b, bSize);
            } else {
                firstAndThirdPhaseCombined<false>(kernel, a, aSize, b, bSize);
            }
        }

    private:
        int threadsNum;
        common::StaggeredStickyBraidMultiplication *staggeredStickyBraidMultiplication;

        /**
          * Computes the kernel of semi-local lcs solution for given two strings with antidiagonal pattern using Open MP where
          * 1st and 3rd phases are merged together so there less syncronizations
        **/
        template<bool Parallel>
        void firstAndThirdPhaseCombined(common::Permutation &matrix, const Input *a, int aSize, const Input *b, int bSize) {
            if (aSize > bSize) {
                auto p = common::Permutation(aSize + bSize, aSize + bSize);
                firstAndThirdPhaseCombined<Parallel>(p, b, bSize, a, aSize);
                matrix = common::Permutation(aSize + bSize, aSize + bSize);
                fillPermutationBa(p, matrix, aSize, bSize);
                return;
            }

            auto m = aSize;
            auto n = bSize;

            auto size = m + n;
            auto *strandMap = new Input[size + 2 * (m - 1)];
            auto thirdPhaseMapSize = m * 2 - 2;
            auto thirdPhaseMap = strandMap + size;

            auto p = common::Permutation(m + n, m + n);
            auto q = common::Permutation(thirdPhaseMapSize, thirdPhaseMapSize);

            auto offset = n - (m - 1);
            auto *aReverse = new Input[m];

#pragma omp parallel num_threads(threadsNum)  default(none) shared(a, aReverse, b, strandMap, size, m, n, matrix, p, q, offset, thirdPhaseMap, thirdPhaseMapSize)
            {

                this->fillAReverse(a, aReverse, m);
                int inThirdPhase = m - 1;
                //    init phase
#pragma omp for simd schedule(static) nowait
                for (int k = 0; k < (m + n); ++k) {
                    strandMap[k] = k;
                }

#pragma omp for simd schedule(static) nowait
                for (int k = 0; k < thirdPhaseMapSize; ++k) {
                    if (k < m - 1) {
                        thirdPhaseMap[k] = 2 * k;
                    } else {
                        thirdPhaseMap[k] = (k - (m - 1)) * 2 + 1;
                    }
                }
#pragma omp barrier


                for (int diag_number = 0; diag_number < m - 1; ++diag_number) {


#pragma omp for simd schedule(static) nowait
                    for (int posInDiag = 0; posInDiag < inThirdPhase; ++posInDiag) {

                        auto top_edge = diag_number + posInDiag;
                        auto left_strand = thirdPhaseMap[posInDiag];
                        auto top_strand = thirdPhaseMap[m - 1 + top_edge];
                        bool r = aReverse[posInDiag] == b[offset + top_edge] || (left_strand > top_strand);

                        if (WithIf) {
                            if (r) std::swap(thirdPhaseMap[posInDiag], thirdPhaseMap[m - 1 + top_edge]);
                        } else {
                            thirdPhaseMap[posInDiag] = (left_strand & (r - 1)) | ((-r) & top_strand);
                            thirdPhaseMap[m - 1 + top_edge] = (top_strand & (r - 1)) | ((-r) & left_strand);
                        }


                    }

#pragma omp for simd schedule(static)
                    for (int posInDiag = inThirdPhase; posInDiag < m; ++posInDiag) {
                        auto topEdge = diag_number + posInDiag + 1 - m;
                        auto leftStrand = strandMap[posInDiag];
                        auto topStrand = strandMap[m + topEdge];
                        bool r = aReverse[posInDiag] == b[topEdge] || (leftStrand > topStrand);

                        if constexpr(WithIf) {
                            if (r) if (r) std::swap(strandMap[posInDiag], strandMap[m + topEdge]);
                        } else {
                            strandMap[posInDiag] = (leftStrand & (r - 1)) | ((-r) & topStrand);
                            strandMap[m + topEdge] = (topStrand & (r - 1)) | ((-r) & leftStrand);
                        }


                    }
                    inThirdPhase--;
                }

                //phase 2
                auto topEdge = m;
                for (int j = 0; j < offset; ++j) {
                    this->template antiDiagonalComputation<Parallel>(strandMap, aReverse, b, m, 0, topEdge, 0, j);
                    topEdge++;
                }

#pragma omp for simd schedule(static) nowait
                for (int l = 0; l < m; l++) {
                    if (l == m - 1) {
                        p.set(strandMap[l], n + l);
                    } else {
                        p.set(strandMap[l], l * 2 + offset);
                    }
                }


#pragma omp for simd schedule(static) nowait
                for (int r = m; r < m + n; r++) {
                    if ((r - m) < offset) {
                        p.set(strandMap[r], r - m);
                    } else {
                        p.set(strandMap[r], (r - m - offset + 1) * 2 + offset - 1);
                    }
                }

#pragma omp for simd schedule(static) nowait
                for (int l = 0; l < m - 1; l++) {
                    q.set(thirdPhaseMap[l], m - 1 + l);
                }


#pragma omp for simd schedule(static) nowait
                for (int r = m - 1; r < m + m - 2; r++) {
                    q.set(thirdPhaseMap[r], r - (m - 1));
                }

#pragma omp barrier
            }

            this->staggeredStickyBraidMultiplication->glueingPartToWhole(p, q, offset, 1, matrix);

            delete[] strandMap;
            delete[] aReverse;
        }

        /**
        *  see theorem 5.21
        * Allows get P_{b,a} when you have P_{a,b}
        */
        void fillPermutationBa(const common::Permutation &ab, common::Permutation &ba, int m, int n) {
            ba.resetAll();
            for (int i = 0; i < ab.rows; ++i) {
                auto col = ab.getColByRow(i);
                if (col != NOPOINT) ba.set(n + m - 1 - i, m + n - 1 - col);
            }
        }


    };

    template<class Input, bool WithIf, bool UseSumBound, bool UseDepth>
    class OpenMPHybridRecursive : public SemiLocalStrategy<Input> {
    public:
        OpenMPHybridRecursive(int parallelThreshold, common::StaggeredStickyBraidMultiplication *staggeredStickyBraidMultiplicator, int depthThreshold,
                              int sumThreshold, AntidiagonalIterativeCombing<Input, WithIf> *iterativeCombingStrat) :
                SemiLocalStrategy<Input>(), staggeredStickyBraidMultiplication(staggeredStickyBraidMultiplicator),
                sumBound(sumThreshold), depth(depthThreshold), parallelDepth(parallelThreshold), iterativeCombing(iterativeCombingStrat) {}


        void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept override {
            hybrid(kernel, a, aSize, b, bSize, parallelDepth, depth, sumBound);
        }

    private:
        int parallelDepth;
        int sumBound;
        int depth;
        common::StaggeredStickyBraidMultiplication *staggeredStickyBraidMultiplication;
        AntidiagonalIterativeCombing<Input, WithIf> *iterativeCombing;

        /**
         * Hybrid appoarch of recursive and iterative combing. Firstly, follows the recursive structure, then switches to the iterative combing
         */
        void hybrid(common::Permutation &perm, const Input *a, int m, const Input *b, int n, int parallelThreshold, int depthThreshold, int sumThreshold) {

            if constexpr(UseSumBound) {
                if (m + n <= sumThreshold) {
                    iterativeCombing->compute(a, m, b, n, perm);
                    return;
                }
            }
            if constexpr(UseDepth) {
                if (depthThreshold <= 0) {
                    iterativeCombing->compute(a, m, b, n, perm);
                    return;
                }
            }
            if constexpr(!UseDepth && !UseSumBound) {
                if (n == 1 && m == 1) {
                    perm = common::Permutation(2, 2);
                    if (a[0] == b[0]) {
                        perm.set(0, 0);
                        perm.set(1, 1);
                    } else {
                        perm.set(1, 0);
                        perm.set(0, 1);
                    }
                    return;
                }
            }

            if (n > m) {
                auto n1 = n / 2;
                auto b1 = b;
                auto b2 = b + n1;

                auto subtreeL = common::Permutation(n1 + m, n1 + m);
                auto subtreeR = common::Permutation(n - n1 + m, n - n1 + m);

                if (parallelThreshold > 0) {
#pragma omp parallel num_threads(2)
                    {
#pragma omp single nowait
                        {
#pragma omp task
                            hybrid(subtreeL, a, m, b1, n1, parallelThreshold - 1, depthThreshold - 1, sumThreshold);

#pragma omp task
                            hybrid(subtreeR, a, m, b2, n - n1, parallelThreshold - 1, depthThreshold - 1, sumThreshold);
                        }

                    }
#pragma omp taskwait
                } else {
                    hybrid(subtreeL, a, m, b1, n1, parallelThreshold - 1, depthThreshold - 1, sumThreshold);
                    hybrid(subtreeR, a, m, b2, n - n1, parallelThreshold - 1, depthThreshold - 1, sumThreshold);
                }
                staggeredStickyBraidMultiplication->staggeredStickyMultiplication<false>(subtreeL, subtreeR, m, perm);
            } else {

                auto m1 = m / 2;
                auto a1 = a;
                auto a2 = a + m1;

                auto subtreeL = common::Permutation(m1 + n, m1 + n);
                auto subtreeR = common::Permutation(m - m1 + n, m - m1 + n);


                if (parallelThreshold > 0) {
#pragma omp parallel num_threads(2)
                    {
#pragma omp single nowait
                        {
#pragma omp task
                            hybrid(subtreeL, a1, m1, b, n, parallelThreshold - 1, depthThreshold - 1, sumThreshold);
#pragma omp task
                            hybrid(subtreeR, a2, m - m1, b, n, parallelThreshold - 1, depthThreshold - 1, sumThreshold);
                        }
                    }
#pragma omp taskwait
                } else {
                    hybrid(subtreeL, a1, m1, b, n, parallelThreshold - 1, depthThreshold - 1, sumThreshold);
                    hybrid(subtreeR, a2, m - m1, b, n, parallelThreshold - 1, depthThreshold - 1, sumThreshold);
                }
                staggeredStickyBraidMultiplication->staggeredStickyMultiplication<true>(subtreeL, subtreeR, n, perm);
            }
        }
    };

    /**
     * To hvae boost use uint
     * @tparam Input
     * @tparam WithIf
     * @tparam UseSumBound
     * @tparam UseDepth
     */
    template<class Input, bool WithIf>
    class OpenMPHybridIterative : public SemiLocalStrategy<Input> {
    public:
        OpenMPHybridIterative(int numThreads, int colsPerBlock, int rowsPerBlock, common::StaggeredStickyBraidMultiplication *staggeredStickyBraidMultiplicator,
                              AntidiagonalIterativeCombing<Input, WithIf> *iterativeCombingStrat) :
                SemiLocalStrategy<Input>(), staggeredStickyBraidMultiplication(staggeredStickyBraidMultiplicator), threads(numThreads),
                iterativeCombing(iterativeCombingStrat), maxColsInBlock(colsPerBlock), maxRowsInBlock(rowsPerBlock) {}

        void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept override {
            /**
             * Heuristic for spliting is as follows.
             * Set as small as possible blocks in one dimenstion bounded by either 32000 or m.
             * If m is greater then  small_m is  m / 32k
             *
             * For other dimenstion. We split equally work among threads. So they get same task with same length
             */

            int colsPerBlock = maxColsInBlock;
            int rowsInBlock = maxRowsInBlock;
            int totalInBlock = colsPerBlock + rowsInBlock;

            int mSmall;
            int nSmall;

            // if m is less then bound then n_small would be 1
            if (rowsInBlock >= aSize) {
                rowsInBlock = aSize;
                mSmall = 1;
            } else {
                // else we check how many blocks in row we need with size rows_per_block
                // and adjust rows_per_block to equal partion
                int nums = int(ceil((1.0 * aSize) / rowsInBlock));
                rowsInBlock = int(ceil(1.0 * aSize / nums));
                mSmall = ceil((aSize * 1.0) / rowsInBlock);
            }


            // if fits then equally split between threads
            if (colsPerBlock >= bSize) {
                colsPerBlock = int(ceil((bSize * 1.0) / threads));
                nSmall = threads;
            } else {
                int cols_per_thread = int(ceil((1.0 * bSize) / threads));
                int rest = (totalInBlock - rowsInBlock);

                if (cols_per_thread < rest) {
                    nSmall = threads;
                } else {
                    nSmall = threads * ceil(1.0 * cols_per_thread / rest);
                }
            }

            hybridIterative(kernel, a, aSize, b, bSize, mSmall, nSmall, threads);
        }

    private:
        int threads;
        int maxColsInBlock;
        int maxRowsInBlock;
        common::StaggeredStickyBraidMultiplication *staggeredStickyBraidMultiplication;
        AntidiagonalIterativeCombing<Input, WithIf> *iterativeCombing;


        /**
         * The hybrid appaorch with down to top apparoch. No recustion. Fixed amount of iteartive combing problems that further are merged via
         * sticky braid multiplication.
         **/
        void hybridIterative(common::Permutation &perm, const Input *a, int m, const Input *b, int n, int smallM, int smallN, int threadsNum = 1) {

            int colsPerBlock = ceil(double(n) / smallN);
            int rowsPerBlock = ceil(double(m) / smallM);
            smallN = std::ceil(double(n) / colsPerBlock);
            smallM = std::ceil(double(m) / rowsPerBlock);


            int numTasks = smallM * smallN;

            auto tasks = new common::Permutation[numTasks];
            auto tasksNextIter = new common::Permutation[numTasks];

#pragma omp parallel  master taskloop num_threads(threadsNum)
            for (int i = 0; i < numTasks; i++) {
                int startCol = (i % smallN) * colsPerBlock;
                int endCol = std::min(startCol + colsPerBlock, n);

                int startRow = (i / smallN) * rowsPerBlock;
                int endRow = std::min(startRow + rowsPerBlock, m);

                // the edge blocks may need do extra work
                if ((i % smallN) == smallN - 1) endCol = n;


                if ((i / smallN) == smallM - 1) endRow = m;


                int sizeBlockB = endCol - startCol;
                int sizeBlockA = endRow - startRow;

                auto b_loop = b + startCol;
                auto a_loop = a + startRow;

                auto matrix = common::Permutation(sizeBlockA + sizeBlockB, sizeBlockA + sizeBlockB);

                iterativeCombing->compute(a_loop, sizeBlockA, b_loop, sizeBlockB, matrix);
                matrix.m = sizeBlockA;
                matrix.n = sizeBlockB;

                tasks[i] = matrix;
            }


            auto nextJobs = tasks;
            auto currentJobs = tasksNextIter;

            int steps = ceil(log2(smallN)) + ceil(log2(smallM));


            int blockRows = rowsPerBlock;
            int blockCols = colsPerBlock;


            double total = 0;
            for (int j = 0; j < steps; ++j) {


                bool isReductionInRow = smallN > smallM;

                //TODO: heuristic now we choose to glue by the large one (seems to work) best. Specifically reduction in a row if on next step row is still big
                //TODO: following is working worse but it copipies logic of recursion: is_reduction_in_row = 2 * block_cols >= block_rows
                if (smallN > 1 && smallM > 1) isReductionInRow = blockRows >= 2 * blockCols;


                if (isReductionInRow) {
                    blockCols *= 2;
                } else {
                    blockRows *= 2;
                }


                auto newCols = isReductionInRow ? int(ceil(smallN / 2.0)) : smallN;
                auto newRows = !isReductionInRow ? int(ceil(smallM / 2.0)) : smallM;

                auto tmp = currentJobs;
                currentJobs = nextJobs;

                nextJobs = tmp;


#pragma omp parallel master taskloop num_threads(threadsNum)
                for (int i = 0; i < newCols * newRows; i++) {

                    auto curRow = i / newCols;
                    auto curCol = i % newCols;


                    common::Permutation p;
                    common::Permutation q;

                    if (isReductionInRow) {
                        p = currentJobs[curRow * smallN + 2 * curCol];

                        if ((2 * curCol + 1) >= smallN) {
                            nextJobs[i] = p;
                        } else {
                            q = currentJobs[curRow * smallN + 2 * curCol + 1];
                            auto product = common::Permutation(p.m + q.n + p.n, p.m + q.n + p.n);

                            staggeredStickyBraidMultiplication->template staggeredStickyMultiplication<false>(p, q, p.m, product);
                            product.m = p.m;
                            product.n = p.n + q.n;

                            nextJobs[i] = product;
                        }

                    } else {
                        p = currentJobs[2 * curRow * smallN + curCol];

                        if ((2 * curRow + 1) >= smallM) {
                            nextJobs[i] = p;
                        } else {
                            q = currentJobs[(2 * curRow + 1) * smallN + curCol];
                            auto product = common::Permutation(p.m + q.m + p.n, p.m + q.m + p.n);
                            staggeredStickyBraidMultiplication->template staggeredStickyMultiplication<true>(p, q, p.n, product);
                            product.m = p.m + q.m;
                            product.n = p.n;
                            nextJobs[i] = product;

                        }
                    }

                }

                smallN = newCols;
                smallM = newRows;

            }

            perm = std::move(nextJobs[0]);
        }

    };

}