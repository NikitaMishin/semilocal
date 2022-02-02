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
    SimpleIterativeCombing<int, true> solver2;
    compareTestEqualStrings(solver1, solver2, start, end, step, high);
}

TEST_F(SemiLocalLCSTest, NaiveVsOpenMPCorrectness) {
    int start = 4900;
    int end = 5000;
    int step = 11;
    int high = 11;
    SimpleIterativeCombing<int, false> solver1;

    for (int i = 1; i <  5; ++i) {
        OpenMPIterativeCombing<int,false> solver2(i);
        compareTestEqualStrings(solver1, solver2, start, end, step, high);
    }
 }

TEST_F(SemiLocalLCSTest, OpenMPSpeedIf) {
    int start = 100000;
    int end = 100001;
    int step = 11;
    int high = 150;
    OpenMPIterativeCombing<int, true> solver1(1);
    compareTestEqualStrings(solver1, solver1, start, end, step, high);
}

TEST_F(SemiLocalLCSTest, OpenMpSpeedNoIf) {
    int start = 100000;
    int end = 100001;
    int step = 11;
    int high = 150;
    OpenMPIterativeCombing<int, false> solver1(1);
    compareTestEqualStrings(solver1, solver1, start, end, step, high);
}