#pragma once


#include "../../common/include/Permutations/Matrices.h"
#include "../SemiLocal/SemiLocalLCS.h"


namespace lcs {
    template<class Input>
    class LCSStrategy {
    public:
        virtual int compute(const Input *a, int sizeA, const Input *b, int sizeB);
    };

    template<class Input>
    class DummyLCS : public LCSStrategy<Input> {
    public:
        int compute(const Input *a, int sizeA, const Input *b, int sizeB) override {
            int arr[sizeA + 1][sizeB + 1];
            auto m = sizeA + 1;
            auto n = sizeB + 1;

            for (auto j = 0; j < n; j++) arr[0][j] = 0;
            for (auto i = 0; i < m; i++) arr[i][0] = 0;

            for (int i = 1; i < m; ++i) {
                for (int j = 1; j < n; ++j) {
                    arr[i][j] = std::max(std::max(arr[i - 1][j], arr[i][j - 1]), (a[i - 1] == b[j - 1]) ? arr[i - 1][j - 1] + 1 : arr[i - 1][j - 1]);
                }
            }
            return arr[m - 1][n - 1];
        }
    };

    template<class Input>
    class SequentialMemLCS : public LCSStrategy<Input> {
    public:
        int compute(const Input *a, int sizeA, const Input *b, int sizeB) override {

            const int *inputA;
            const int *inputB;
            int m, n;

            if (a > b) {
                m = sizeA + 1;
                n = sizeB + 1;
                inputA = a;
                inputB = b;
            } else {
                n = sizeA + 1;
                m = sizeB + 1;
                inputB = a;
                inputA = b;
            }

            auto prevRow = new int[n];
            auto curRow = new int[n];
            for (int i = 0; i < n; ++i) {
                curRow[i] = 0;
                prevRow[i] = 0;
            }

            for (int i = 1; i < m; ++i) {
                auto l = 0;
                for (int j = 1; j < n; ++j) {
                    curRow[j] = std::max(
                            std::max(prevRow[j], l),
                            (inputA[i - 1] == inputB[j - 1]) ? prevRow[j - 1] + 1 : prevRow[j - 1]
                    );
                    l = curRow[j];

                }
                std::swap(prevRow, curRow);
            }
            auto res = prevRow[n - 1];
            delete[] prevRow;
            delete[] curRow;
            return res;
        }
    };

    template<class Input>
    class OpenMPSIMDAntidiagonalLCS : public LCSStrategy<Input> {
    public:
        //TODO check performance with classes &  simplify code
        int compute(const Input *a, int sizeA, const Input *b, int sizeB) override {
            if (sizeA > sizeB) {
                return compute(b, sizeB, a, sizeA);
            }

            if (sizeA == 1 && sizeB == 1) {
                return a[0] == b[0] ? 1 : 0;
            }

            auto diagonalSize = 1 + std::min(sizeA, sizeB);
            auto a1 = new int[diagonalSize];
            auto a2 = new int[diagonalSize];
            auto a3 = new int[diagonalSize];

            auto startI = 0;
            auto startJ = 0;
            auto min = std::min(sizeA, sizeB);
            auto numDiag = sizeA + sizeB;
            auto totalSameLengthDiag = numDiag - (min + 1) - min;

#pragma omp simd
            for (int k = 0; k < diagonalSize; ++k) {
                a3[k] = 0;
                a2[k] = 0;
                a1[k] = 0;
            }


            startI--;
            // fill upper square
            for (int k = 2; k <= diagonalSize; ++k, startI++) {
                a3[0] = 0;
                a3[k - 1] = 0;
#pragma omp simd
                for (int i = 1; i < k - 1; ++i) {
                    a3[i] = std::max(
                            std::max(a2[i], a2[i - 1]),
                            (a[startI + 1 - i] == b[-1 + i]) ? 1 + a1[i - 1] : a1[i - 1]
                    );
                }
                std::swap(a1, a2);
                std::swap(a3, a2);
            }

            // phase 2:: fill
            if (sizeA >= sizeB) {
                //        same pattern
                for (int k = 0; k < totalSameLengthDiag; ++k, startI++) {
                    a3[0] = 0;

#pragma omp simd
                    for (int i = 1; i < diagonalSize; ++i) {
                        a3[i] = std::max(
                                std::max(a2[i], a2[i - 1]),
                                (a[startI + 1 - i] == b[i - 1]) ? 1 + a1[i - 1] : a1[i - 1]
                        );
                    }
                    std::swap(a1, a2);
                    std::swap(a3, a2);
                }

            }

            //      special case when:
            //      a==b => |a1| = c-1 , |a2| = c, |a3|= c-1  or
            //      a>b  => |a1| = c, |a2| = c, |a3| = c-1
            //      a<b ->  |a1| = c - 1, |a2| = c, |a3| = c

            if (sizeA < sizeB) {
                a3[diagonalSize - 1] = 0;
            }

#pragma omp simd
            for (int i = 0; i < diagonalSize - 1; ++i) {
                a3[i] = std::max(
                        std::max(a2[i], a2[i + 1]),
                        (a[startI - i] == b[i]) ? 1 + a1[i] : a1[i]
                );
            }
            startJ++;
            std::swap(a1, a2);
            std::swap(a3, a2);

            if (sizeA < sizeB) {
//        since special case then -1
                for (int k = 0; k < totalSameLengthDiag; ++k, startJ++) {

                    a3[diagonalSize - 1] = 0;
#pragma omp simd
                    for (int i = 0; i < diagonalSize - 1; ++i) {
                        a3[i] = std::max(
                                std::max(a2[i], a2[i + 1]),
                                (a[startI - i] == b[startJ + i]) ? 1 + a1[i + 1] : a1[i + 1]
                        );
                    }
                    std::swap(a1, a2);
                    std::swap(a3, a2);
                }
            }

            if (sizeA >= sizeB) diagonalSize -= 1;


//    phase 3
//    pattern a3[i] = max(a1[i+1],a2[i],a2[i-1])
            for (int size = diagonalSize - 1; size > 1; size--, startJ++) {

                for (int i = 0; i < size; ++i) {
                    a3[i] = std::max(
                            std::max(a2[i], a2[i + 1]),
                            (a[startI - i] == b[startJ + i]) ? 1 + a1[i + 1] : a1[i + 1]
                    );

                }

                std::swap(a1, a2);
                std::swap(a3, a2);
            }


            auto res = std::max(std::max(a2[0], a2[1]), (a[sizeA - 1]) == b[sizeB - 1] ? 1 + a1[1] : a1[1]);
            delete[]a1;
            delete[]a2;
            delete[]a3;
            //  need to calculate last one cell
            return res;
        }

    };

    template<class Input>
    class LCSBySemiLocal : public LCSStrategy<Input> {

    public:

        explicit LCSBySemiLocal(semi_local::SemiLocalStrategy<Input> *semiLocalSolver) : solver(semiLocalSolver) {}

        int compute(const Input *a, int sizeA, const Input *b, int sizeB) override {
            common::Permutation kernel;
            solver->compute(a, sizeA, b, sizeB, kernel);

            auto acc = 0;
            for (int startOfStrand = 0; startOfStrand < sizeA; ++startOfStrand) {
                auto endOfStrand = kernel.getColByRow(startOfStrand);
                if (sizeB <= endOfStrand) acc++;
            }
            return sizeA - acc;
        }

    private:
        semi_local::SemiLocalStrategy<Input> *solver;

    };


}