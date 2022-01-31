
#include <gtest/gtest.h>
#include "../../common/include/Permutations/Matrices.h"
#include "../../common/include/Permutations/MongeMultuplicationStrategy.h"

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
};

TEST_F(SteadyAntTest, A_LOT_OF_GENERATED_TESTS) {
    using namespace common;

    NaiveBraidMultiplication solver;

    PrecalcMap map;
    PrecalcMap emptyMap;
    NaiveBraidMultiplication tmpSolver;
    StickyBraidMultiplication::buildPrecalcMap(&tmpSolver, map, 5);

    for (int size = 100; size < 300; ++size) {
        Permutation P(size, size);
        Permutation Q(size, size);
        Permutation expected(size, size);
        Permutation actual(size, size);
        fillPermutationMatrix(P, size, size);
        fillPermutationMatrix(Q, size, size);

        solver.multiply(P, Q, expected);
        SimpleStickyBraidMultiplication stickySolver(map);

        stickySolver.multiply(P, Q, actual);


        ASSERT_EQ(expected, actual);


    }
}


TEST_F(SteadyAntTest, MultiplicatonCorrecntess) {
    using namespace common;
    Permutation P(smallSize, smallSize);
    Permutation Q(smallSize, smallSize);
    Permutation R(smallSize, smallSize);
    fillPermutationMatrix(P, smallSize, smallSize);
    fillPermutationMatrix(Q, smallSize, smallSize);
    NaiveBraidMultiplication solver;
    solver.multiply(P, Q, R);

    PrecalcMap emptyMap;
    SimpleStickyBraidMultiplication stickySolver(emptyMap);
    Permutation R2(smallSize, smallSize);
    stickySolver.multiply(P, Q, R2);

    ASSERT_EQ(R2, R);

}