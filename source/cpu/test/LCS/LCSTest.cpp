#include <gtest/gtest.h>
#include <omp.h>
#include <chrono>
#include "../../include/BitwiseLCS/BitParalellLCS.h"
#include "../../common/include/Permutations/Matrices.h"
#include "../../common/include/Permutations/MongeMultuplicationStrategy.h"
#include "../../include/SemiLocal/SemiLocalLCS.h"
#include "../../include/LCS/LCS.h"

using namespace semi_local;
using namespace lcs;

class LCSTest : public ::testing::Test {
public:

    LCSTest() {
        srand(seed);
    }

    template<class Input>
    void fillStringRandomly(Input *string, int high, int sz) {
        for (int i = 0; i < sz; ++i) string[i] = abs(rand()) % high;
    }

protected:
    template<typename T>
    struct Tag {
        typedef T type;
    };

    int seed = 0;


    template<class Input>
    void compareTestEqualStrings(LCSStrategy<Input> &solver1, LCSStrategy<Input> &solver2, int start, int end, int step, int high) {
        using namespace common;
        for (int sizeA = start; sizeA < end; sizeA += step) {
            for (int sizeB = start; sizeB < end; sizeB += step) {
                Permutation expected;
                Permutation actual;
                auto A = new Input[sizeA];
                auto B = new Input[sizeB];
                fillStringRandomly(A, high, sizeA);
                fillStringRandomly(B, high, sizeB);
                auto res1 = solver1.compute(A, sizeA, B, sizeB);
                auto res2 = solver2.compute(A, sizeA, B, sizeB);
                if (res1 != res2) {
                    std::cout << sizeA << "," << sizeB << std::endl;
                    for (int i = 0; i < sizeA; ++i) std::cout << A[i];
                    std::cout << std::endl;
                    for (int i = 0; i < sizeB; ++i) std::cout << B[i];
                    ASSERT_EQ(res1, res2);
                }
                delete[] A;
                delete[] B;
            }
        }
    }


    template<typename T,bool useNaive, bool formula1>
    void testOur(int start_, int end_, int alphabetSize, int thds = 4,bool justPerformance = false) {
        auto bits = 8 * sizeof(T);
        int start = start_ * bits;
        int end = end_ * bits + 1;
        int step = bits;
        int high = alphabetSize;
        SequentialMemLCS<int> solver1{};
        auto getSolver = [&](){
            if constexpr(useNaive) {
                return bit_parallel::lcs::BitSemiLocalSticky<T, int, true, true, true>(thds);
            } else if constexpr (formula1) {
                return bit_parallel::lcs::BitSemiLocalSticky<T, int, true, true, false, false>(thds);
            } else {
                return bit_parallel::lcs::BitSemiLocalSticky<T, int, true, true, false, true>(thds);
            }
        };
        auto solver2 = getSolver();
        if (justPerformance) {
                compareTestEqualStrings(solver2, solver2, start, end, step, high);
        } else {
                compareTestEqualStrings(solver1, solver2, start, end, step, high);
        }
    };


};


TEST_F(LCSTest, DummyVsNaiveBySemi) {
    int start = 5;
    int end = 500;
    int step = 1;
    int high = 15;
    DummyLCS<int> solver1{};
    SimpleIterativeCombing<int, false> solverSemi;
    LCSBySemiLocal<int> solver2(&solverSemi);

    compareTestEqualStrings(solver1, solver2, start, end, step, high);
}

TEST_F(LCSTest, DummyVsMemDummy) {
    int start = 5;
    int end = 500;
    int step = 1;
    int high = 15;
    DummyLCS<int> solver1{};
    SequentialMemLCS<int> solver2{};

    compareTestEqualStrings(solver1, solver2, start, end, step, high);
}

TEST_F(LCSTest, DummyVsAntidiagonal) {
    int start = 5;
    int end = 500;
    int step = 1;
    int high = 50;
    DummyLCS<int> solver1{};
    OpenMPSIMDAntidiagonalLCS<int> solver2{};

    compareTestEqualStrings(solver1, solver2, start, end, step, high);
}


TEST_F(LCSTest, MemDummyVsHyyroSimple) {
    int start = 16;
    int end = 600;
    int step = 1;
    int high = 50;
    SequentialMemLCS<int> solver1{};
    bit_parallel::lcs::Hyyro<uint32_t, int> solver2{};

    compareTestEqualStrings(solver1, solver2, start, end, step, high);
}

TEST_F(LCSTest, MemDummyVsBitwiseBinary) {
    auto ts1 = std::chrono::high_resolution_clock::now();
    testOur<uint8_t,false, false>(1, 32, 2);
    testOur<uint16_t,false, false>(1, 32, 2);
    testOur<uint32_t,false, false>(1, 32, 2);
    testOur<uint64_t,false, false>(1, 32, 2);
    auto ts2 = std::chrono::high_resolution_clock::now();
    testOur<uint8_t,true, true>(1, 32, 2);
    testOur<uint16_t,true, true>(1, 32, 2);
    testOur<uint32_t,true, true>(1, 32, 2);
    testOur<uint64_t,true, true>(1, 32, 2);
    auto ts3 = std::chrono::high_resolution_clock::now();
    testOur<uint8_t,true, false>(1, 32, 2);
    testOur<uint16_t,true, false>(1, 32, 2);
    testOur<uint32_t,true, false>(1, 32, 2);
    testOur<uint64_t,true, false>(1, 32, 2);
    auto ts4 = std::chrono::high_resolution_clock::now();
    std::cout << (ts2 - ts1).count() << std::endl;
    std::cout << (ts3 - ts2).count() << std::endl;
    std::cout << (ts4 - ts3).count() << std::endl;
}


TEST_F(LCSTest, PerfomanceBitwiseOurs) {
    auto sz = 32 * 10 * 10*10;
    auto ts1 = std::chrono::high_resolution_clock::now();
    testOur<uint32_t,true, false>(sz, sz, 2,  4, true);
    auto ts2 = std::chrono::high_resolution_clock::now();
    testOur<uint32_t,false, true>(sz, sz, 2,  4, true);
    auto ts3 = std::chrono::high_resolution_clock::now();
    testOur<uint32_t,false, false>(sz, sz, 2,  4, true);
    auto ts4 = std::chrono::high_resolution_clock::now();
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(ts2 - ts1).count() / 2 << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(ts3 - ts2).count() / 2 << std::endl;
    std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(ts4 - ts3).count() / 2 << std::endl;
}



