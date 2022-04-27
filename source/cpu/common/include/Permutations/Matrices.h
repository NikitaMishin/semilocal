#pragma once


#include <memory>

#define NOPOINT (-1) /*NOPOINT indicates that no non-zero element in matrix in specified position*/

namespace common {

    /**
    * Implementation of a permutation matrix is based on two arrays:
    * the  first array is mapping of non-zero entries in rows to its position in cols
    * the second array is mapping of non-zero entries in cols to its position in rows
    * Memory freed outised of class.
    */
    class Permutation  {
    public:
        Permutation() = default;
        ~Permutation() = default;
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
        std::vector<int> arr;
    public:

        Matrix(int row_size, int col_size) {
            rows = row_size;
            cols = col_size;
            arr.resize(row_size * col_size, 0);
        }

        int getElementAt(int row, int col) const { return arr[row * cols + col]; }

        inline void setElementAt(int row, int col, int value) { arr[row * cols + col] = value; }

        void getCrossDiffPermutation(Permutation &perm){
            perm = Permutation(rows - 1, cols - 1);
            for (int i = 0; i < rows - 1; ++i) {
                for (int j = 0; j < cols - 1; ++j) {
                    auto crossDiff = (getElementAt(i, j + 1) + getElementAt(i + 1, j)) - (getElementAt(i, j) + getElementAt(i + 1, j + 1));
                    if (crossDiff == 1) perm.set(i, j);
                }
            }
        }


        int rows;
        int cols;
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



