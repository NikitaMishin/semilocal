#pragma once

#include "../../common/include/Permutations/MongeMultuplicationStrategy.h"
#include "../LCS/LCS.h"

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

        void compute(const Input *a, int aSize, const Input *b, int bSize, common::Permutation &kernel) noexcept override  {
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
            if(threadsNum > 1) {
                stickyBraidMpi<true>(kernel, a, aSize, b, bSize);
            } else {
                stickyBraidMpi<false>(kernel,a,aSize,b,bSize);
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
            matrix = common::Permutation(m+n,m+n);

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
                    this-> template antiDiagonalComputation<Parallel>(strand_map, aReverse, b, m, 0, topEdge, 0, j);
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


}