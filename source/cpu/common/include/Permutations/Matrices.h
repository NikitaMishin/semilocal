#pragma once


#include <memory>

#define NOPOINT (-1) /*NOPOINT indicates that no non-zero element in matrix in specified position*/

namespace common {

    /**
    * Implementation of (sub)permutation matrix is based on two arrays:
    * the  first array is mapping of non-zero entries in rows to its position in cols
    * the second array is mapping of non-zero entries in cols to its position in rows
    * Note that memory management is handled inside class
    */
//    class Permutation {
//    public:
//
//        /**
//        * Initializes permutation matrix of size rowXcol  as zero matrix that have no non-zero entries
//        * @param row
//        * @param col
//        */
//        Permutation(int row, int col) {
//            rows = row;
//            cols = col;
//            rowToCol = new int[row];
//            colToRow = new int[col];
//            for (int i = 0; i < rows; ++i) rowToCol[i] = NOPOINT;
//            for (int i = 0; i < cols; ++i) colToRow[i] = NOPOINT;
//        }
//
//        inline void set(int row, int col)  {
//            rowToCol[row] = col;
//            colToRow[col] = row;
//        }
//
//        inline void reset(int row, int col) {
//            if (rowToCol[col] == col) {
//                rowToCol[col] = NOPOINT;
//                colToRow[row] = NOPOINT;
//            }
//        }
//
//        inline void resetAll() {
//            for (int i = 0; i < rows; ++i) rowToCol[i] = NOPOINT;
//            for (int i = 0; i < cols; ++i) colToRow[i] = NOPOINT;
//        }
//
//        inline int getRowByCol(int col) const  { return colToRow[col]; }
//
//        inline int getColByRow(int row) const  { return rowToCol[row]; }
//
//        int getElementAt(int row, int col) const { return (getColByRow(row) == col) ? 1 : 0; }
//
//        bool operator==(const Permutation &other) const {
//            if (other.rows != rows || other.cols != cols) return false;
//            for (int i = 0; i < rows; ++i) if (getColByRow(i) != other.getColByRow(i)) return false;
//            for (int i = 0; i < cols; ++i) if (getRowByCol(i) != other.getRowByCol(i)) return false;
//            return true;
//        }
//
//
//        /**
//        * Fills input vector with position pairs of non-zero entries in the current matrix aka (row_i,col_i)
//        * @param result std::vector<std::pair<int,int>>
//        * @return void
//        */
//        void toPointsOnGrid(std::vector<std::pair<int, int>> &result) const {
//            for (int i = 0; i < cols; ++i) {
//                auto col = getColByRow(i);
//                if (col != NOPOINT) result.emplace_back(i, col);
//            }
//        }
//
//
//        /**
//         * Initializes permutation matrix of size rowXcol  as matrix that has non-zero entries given by vector points
//         * @param row
//         * @param col
//         * @param points
//         */
//        explicit Permutation(int row, int col, std::vector<std::pair<int, int>> &points) : Permutation(row, col) {
//            for (auto &point: points) {
//                rowToCol[point.first] = point.second;
//                colToRow[point.second] = point.first;
//            }
//        }
//
//        ~Permutation() {
//            if(rowToCol != nullptr) delete[] rowToCol;
//            if(colToRow != nullptr) delete[] colToRow;
//        }
//
//        //TODO after refactoring move to sep instance
//        int m = 0;
//        int n = 0;
//
//    private:
//        int* rowToCol;
//        int* colToRow;
//        int rows;
//        int cols;
//    };

    /**
    * Implementation of a permutation matrix is based on two arrays:
    * the  first array is mapping of non-zero entries in rows to its position in cols
    * the second array is mapping of non-zero entries in cols to its position in rows
    * Memory freed outised of class.
    */
    class Permutation  {
    public:
        Permutation() = default;
        Permutation(int row, int col) :rows(row), cols(col) {
            rowToCol.resize(row);
            std::fill(rowToCol.begin(), rowToCol.end(), NOPOINT);
            colToRow.resize(col);
            std::fill(colToRow.begin(),colToRow.end(),NOPOINT);
        };


        bool operator==(const Permutation &other) const {
            if (other.rows != rows || other.cols != cols) return false;
            for (int i = 0; i < rows; ++i) if (getColByRow(i) != other.getColByRow(i)) return false;
            for (int i = 0; i < cols; ++i) if (getRowByCol(i) != other.getRowByCol(i)) return false;
            return true;
        }


        void fromPoints(const std::vector<std::pair<int, int>> &points) {
            rows  = points.size();
            cols = points.size();
            rowToCol.resize(points.size());
            colToRow.resize(points.size());
            for (auto &point: points) {
                rowToCol[point.first] = point.second;
                colToRow[point.second] = point.first;
            }
        }

        void fromHash() {

        }

        /**
        *
        * Fills input vector with position pairs of non-zero entries in the current matrix aka (row_i,col_i)
        * @param result std::vector<std::pair<int,int>>
        * @return void
        * NOTICE for square matrix only
        */
        void toPoints(std::vector<std::pair<int, int>> &result) const {
            for (int i = 0; i < cols; ++i) {
                auto col = getColByRow(i);
                if (col != NOPOINT) result.emplace_back(i, col);
            }
        }

        inline void set(int row, int col)  {
            rowToCol[row] = col;
            colToRow[col] = row;
        }

        inline void reset(int row, int col)  {
            if (rowToCol[col] == col) {
                rowToCol[col] = NOPOINT;
                colToRow[row] = NOPOINT;
            }
        }

        int getElementAt(int row, int col) const { return (getColByRow(row) == col) ? 1 : 0; }


        inline void resetAll()  {
            std::fill(rowToCol.begin(),rowToCol.end(),NOPOINT);
            std::fill(colToRow.begin(),colToRow.end(),NOPOINT);
        }

        inline int getRowByCol(int col) const { return colToRow[col]; }

        inline int getColByRow(int row) const { return rowToCol[row]; }




        //TODO after refactoring move to sep instance
        int m = 0;
        int n = 0;
        int rows;
        int cols;
    private:
        std::vector<int> rowToCol;
        std::vector<int> colToRow;
    };

    /**
    * The class presents a simple 2d matrix with integer value entities. Used primarily for testing
    */
    class Matrix  {
    private:
        int *arr;
        int rows;
        int cols;
    public:

        Matrix(int row_size, int col_size) {
            rows = row_size;
            cols = col_size;
            arr = new int[row_size * col_size];
            for (int i = 0; i < row_size * col_size; ++i) arr[i] = 0;
        }

        int getElementAt(int row, int col) const { return arr[row *cols + col]; }

        inline void set_element_at(int row, int col, int value) { arr[row * cols + col] = value; }


        //

        ~Matrix() {
            if(arr!= nullptr) delete[] arr;
        }
    };

    /**
    * Fill given zero permutation matrix randomly by given seed
    * @param m
    * @param rowSize
    * @param colSize
    * @param seed
    */
    inline void fillPermutationMatrix(Permutation &m, int rowSize, int colSize) {

        auto isUsed = new bool[colSize];
        for (int i = 0; i < colSize; ++i) isUsed[i] = false;

        auto activePts = colSize;

        for (int row = 0; row < rowSize && activePts > 0; ++row) {
            while (true) {
                auto col = abs(rand()) % colSize;
                if (!isUsed[col]) {
                    m.set(row, col);
                    isUsed[col] = true;
                    break;
                }
            }
            activePts--;
        }

        delete[] isUsed;
    }



    template<typename Perm>
    inline void debugPrint(const Perm& perm,std::ostream &os)  {
        for (int i = 0; i < perm.rows; ++i) {
            for (int j = 0; j < perm.cols; ++j) {
                os << perm.getElementAt(i, j);
            }
            os << std::endl;
        }
    }

}


///**
// * Fill given zero permutation matrix randomly by given seed
// * @param m
// * @param row_size
// * @param col_size
// * @param seed
// */
//void fill_permutation_matrix(AbstractPermutation *m, int row_size, int col_size, int seed = 0) {
//
//    auto is_used = new bool[col_size];
//    for (int i = 0; i < col_size; ++i) is_used[i] = false;
//
//    /* initialize random seed: */
//    srand(seed);
//
//    auto active_pts = col_size;
//
//    for (int row = 0; row < row_size && active_pts > 0; ++row) {
//
//        while (true) {
//            auto col = abs(rand()) % col_size;
//            if (!is_used[col]) {
//                m->set(row, col);
//                is_used[col] = true;
//                break;
//            }
//        }
//        active_pts--;
//    }
//
//    delete[] is_used;
//}
//
//
///**
// * Get dominance matrix of specified func operator( could left/right bottom/top arrow) from permutations matrix m
// * @tparam Lambda
// * @param m
// * @param func
// * @param output_dominance_matrix
// */
//template<typename Lambda>
//void get_dominance_matrix(AbstractPermutation &m, Lambda &&func, Matrix *output_dominance_matrix) {
//    auto row_size = m.row_size + 1;
//    auto col_size = m.col_size + 1;
//
//    for (int row = 0; row < row_size; ++row) {
//        for (int col = 0; col < col_size; ++col) {
//            for (int i = 0; i < m.row_size; ++i) {
//                auto row_pos_point = i;
//                auto col_pos_point = m.getColByRow(row_pos_point);
//                if (col_pos_point == NOPOINT) continue;
//                if (func(row_pos_point, col_pos_point, row, col) == true)
//                    output_dominance_matrix->set_element_at(row, col,
//                                                            output_dominance_matrix->getElementAt(row, col) + 1);
//            }
//        }
//    }
//}
//
//
///**
// * Cope elements from permutation matrix from to permutation matrix to
// * @param from
// * @param to
// */
//inline void copy(AbstractPermutation &from, AbstractPermutation &to) {
//    for (int i = 0; i < from.row_size; ++i) {
//        auto col = from.getColByRow(i);
//        to.set(i, col);
//    }
//}


