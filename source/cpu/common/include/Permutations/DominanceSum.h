
#include "Matrices.h"

namespace common {

    enum Arrow {
        TOP_LEFT,
        TOP_RIGHT,
        BOTTOM_LEFT,
        BOTTMO_RIGHT
    };
    enum Move {
        LEFT,
        UP,
        DOWN,
        RIGHT
    };

    /**
    * Provides ability to perform incremental queries
    */
    class IncrementalDominanceSum {


    public:
        /**
         * @tparam arrow
         * @tparam moveAction
         * @param row
         * @param col
         * @param sum
         * @param permMatrix
         * @return
         */
        template<Arrow arrow, Move moveAction>
        static int move(int row, int col, int sum, const Permutation &permMatrix) noexcept {
            if constexpr(arrow == TOP_LEFT) {
                if constexpr(moveAction == LEFT) {
                    if (col == 0) return sum;
                    auto rowCap = permMatrix.getRowByCol(col - 1);
                    if (rowCap >= row && rowCap != NOPOINT) sum++;
                    return sum;
                }
                if constexpr(moveAction == DOWN) {
                    if (row >= permMatrix.rows) return 0;
                    auto colCap = permMatrix.getColByRow(row);
                    if (colCap >= col && colCap != NOPOINT) sum--;
                    return sum;
                }
                if constexpr(moveAction == UP) {
                    if (row == 0) return sum;
                    auto colCap = permMatrix.getColByRow(row - 1);
                    if (colCap >= col && colCap != NOPOINT) sum++;
                    return sum;
                }
                if constexpr(moveAction == RIGHT) {
                    if (col >= permMatrix.cols) return 0;
                    auto rowCap = permMatrix.getRowByCol(col);
                    if (rowCap >= row && rowCap != NOPOINT) sum--;
                    return sum;
                }
            }

            if constexpr(arrow == BOTTMO_RIGHT) {
                if constexpr(moveAction == LEFT) {
                    if (col == 0) return sum;
                    auto rowCap = permMatrix.getRowByCol(col - 1);

                    if (rowCap != NOPOINT && rowCap < row) sum--;
                    return sum;
                }
                if constexpr(moveAction == DOWN) {
                    if (row >= permMatrix.rows) return 0;
                    auto colCap = permMatrix.getColByRow(row);
                    if (colCap < col && colCap != NOPOINT) sum++;
                    return sum;
                }
                if constexpr(moveAction == UP) {
                    if (row == 0) return sum;
                    auto colCap = permMatrix.getColByRow(row - 1);
                    if (colCap < col && colCap != NOPOINT) sum--;
                    return sum;
                }
                if constexpr(moveAction == RIGHT) {
                    if (col >= permMatrix.cols) return 0;
                    auto rowCap = permMatrix.getRowByCol(col);
                    if (rowCap < row && rowCap != NOPOINT) sum++;
                    return sum;
                }
            }

            if constexpr(arrow == TOP_RIGHT) {
                if constexpr(moveAction == LEFT) {
                    if (col == 0) return sum;
                    auto rowCap = permMatrix.getRowByCol(col - 1);
                    if (rowCap >= row && rowCap != NOPOINT) sum--;
                    return sum;
                }
                if constexpr(moveAction == DOWN) {
                    if (row >= permMatrix.rows) return 0;
                    auto colCap = permMatrix.getColByRow(row);
                    if (colCap < col && colCap != NOPOINT) sum--;
                    return sum;
                }
                if constexpr(moveAction == UP) {
                    if (row == 0) return sum;
                    auto colCap = permMatrix.getColByRow(row - 1);
                    if (colCap < col && colCap != NOPOINT) sum++;
                    return sum;
                }
                if constexpr(moveAction == RIGHT) {
                    if (col >= permMatrix.cols) return 0;
                    auto rowCap = permMatrix.getRowByCol(col);

                    if (rowCap >= row && rowCap != NOPOINT) sum++;
                    return sum;
                }
            }
        }
    };

    /**
     * Get dominance matrix of specified func operator( could left/right bottom/top arrow) from permutations matrix m
    * @tparam Lambda
    * @param m
    * @param func
    * @param outputDominanceMatrix
    */
    template<typename Lambda>
    inline void getDominanceSum(const Permutation &m, Lambda &&func, Matrix &outputDominanceMatrix) {
        auto row_size = m.rows + 1;
        auto col_size = m.cols + 1;

        for (int row = 0; row < row_size; ++row) {
            for (int col = 0; col < col_size; ++col) {
                for (int i = 0; i < m.rows; ++i) {
                    auto row_pos_point = i;
                    auto col_pos_point = m.getColByRow(row_pos_point);
                    if (col_pos_point == NOPOINT) continue;
                    if (func(row_pos_point, col_pos_point, row, col) == true)
                        outputDominanceMatrix.set_element_at(row, col, outputDominanceMatrix.getElementAt(row, col) + 1);
                }
            }
        }
    }

};