#include <gtest/gtest.h>
#include <omp.h>
#include "../../common/include/Permutations/Matrices.h"
#include "../../common/include/Permutations/MongeMultuplicationStrategy.h"
#include "../../include/semiLocal/SemiLocalLCS.h"


using namespace semi_local;

class SemiLocalLCSTest : public ::testing::Test {
public:

    SemiLocalLCSTest() {
        srand(seed);
    }

    void fillStringRandomly(int *string, int high, int sz) {
        for (int i = 0; i < sz; ++i) string[i] = abs(rand()) % high;
    }

protected:
    int seed = 0;


    void compareTestEqualStrings(SemiLocalStrategy<int> &solver1, SemiLocalStrategy<int> &solver2, int start, int end, int step, int high) {
        using namespace common;
        for (int sizeA = start; sizeA < end; sizeA += step) {
            for (int sizeB = start; sizeB < end; sizeB += step) {
                Permutation expected;
                Permutation actual;
                auto A = new int[sizeA];
                auto B = new int[sizeB];
                fillStringRandomly(A, high, sizeA);
                fillStringRandomly(B, high, sizeB);
                solver1.compute(A, sizeA, B, sizeB, expected);
                solver2.compute(A, sizeA, B, sizeB, actual);
                ASSERT_EQ(expected, actual);
                delete[] A;
                delete[] B;
            }
        }
    }
};

TEST_F(SemiLocalLCSTest, DummyVsNaiveCorrectness) {
    int start = 5;
    int end = 50;
    int step = 1;
    int high = 5;
    DummyPrefixLCSStrategy<int> solver1(-2);
    SimpleIterativeCombingStrategy<int, true> solver2;
    SimpleIterativeCombingStrategy<int, true> solver3;

    compareTestEqualStrings(solver1, solver2, start, end, step, high);

}
