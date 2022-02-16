#include <gtest/gtest.h>
#include <omp.h>
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
                ASSERT_EQ(solver1.compute(A, sizeA, B, sizeB), solver2.compute(A, sizeA, B, sizeB));
                delete[] A;
                delete[] B;
            }
        }
    }
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
    bit_parallel::lcs::Hyyro<uint32_t ,int>solver2{};

    compareTestEqualStrings(solver1, solver2, start, end, step, high);
}
