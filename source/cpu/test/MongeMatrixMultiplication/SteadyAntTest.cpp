
#include <gtest/gtest.h>
#include <omp.h>
#include "../../common/include/Permutations/Matrices.h"
#include "../../common/include/Permutations/MongeMultuplicationStrategy.h"
#include "../../include/semiLocal/SemiLocalLCS.h"
class SteadyAntTest : public ::testing::Test {
public:

    SteadyAntTest() {
        /* initialize random seed: */
        srand(seed);
    }

    void fillPermMatrixRandomly(common::Permutation &perm, int sz) {
        common::fillPermutationMatrix(perm, sz, sz);
    }

protected:
    common::PrecalcMap emptyMap;
    int smallSize = 7;
    int mediumSize = 10000;
    int largeSize = 1000000;
    int seed = 0;


    void compareTest(common::BraidMultiplicationStrategy & solver1, common::BraidMultiplicationStrategy & solver2,int start,int end,int step) {
        using namespace common;
        for (int size = start; size < end;size+=step) {
            Permutation P(size, size);
            Permutation Q(size, size);
            Permutation expected;
            Permutation actual;
            fillPermutationMatrix(P, size, size);
            fillPermutationMatrix(Q, size, size);

            solver1.multiply(P, Q, expected);
            solver2.multiply(P, Q, actual);
            ASSERT_EQ(expected, actual);
        }
    }
};

TEST_F(SteadyAntTest, NaiveVsSimpleSticky) {
    using namespace common;
    NaiveBraidMultiplication solver;
    PrecalcMap emptyMap;
    SimpleStickyBraidMultiplication tmpSolver(emptyMap);

    for(int i =0;i<5;i++) {
        PrecalcMap map;
        StickyBraidMultiplication::buildPrecalcMap(&tmpSolver, map, i);

        NaiveBraidMultiplication solver1;
        SimpleStickyBraidMultiplication solver2(map);
        compareTest(solver1, solver2, 100, 300, 1);
    }
}


TEST_F(SteadyAntTest, SimpleStickyVsMemory) {
    using namespace common;
    PrecalcMap map;
    NaiveBraidMultiplication tmpSolver;
    StickyBraidMultiplication::buildPrecalcMap(&tmpSolver, map, 5);

    SimpleStickyBraidMultiplication solver1(map);
    SequentialMemoryOptimizedStickBraidMultiplication solver2(map);
    compareTest(solver1,solver2,1000,30000,111);
}

TEST_F(SteadyAntTest, SimpleStickyVsMemoryLong) {
    using namespace common;
    PrecalcMap map;
    NaiveBraidMultiplication tmpSolver;
    StickyBraidMultiplication::buildPrecalcMap(&tmpSolver, map, 5);

    SimpleStickyBraidMultiplication solver1(map);
    SequentialMemoryOptimizedStickBraidMultiplication solver2(map);
    compareTest(solver1,solver2,1000000,1000001,100);
}

TEST_F(SteadyAntTest, StickyOptMemVsOpenMPLong) {
    using namespace common;
    PrecalcMap map;
    NaiveBraidMultiplication tmpSolver;
    StickyBraidMultiplication::buildPrecalcMap(&tmpSolver, map, 5);

    SequentialMemoryOptimizedStickBraidMultiplication solver1(map);
    OpenMPStickyBraid solver2(2,map);
    omp_set_nested(true);
    compareTest(solver1,solver2,1000000,1000001,100);
}

