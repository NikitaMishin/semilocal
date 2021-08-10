//
// Created by garrancha on 10.08.2021.
//

#ifndef CPU_SMAWK_H
#define CPU_SMAWK_H

#include <functional>
#include <stack>


// todo add linear
class RowMinima {
public:
    RowMinima(std::function<double(int, int)> accessor) : _accessor(accessor) {};

    void rowMinima(int m, int n, std::vector<int> &output) {
        for (int i = 0; i < m; ++i) output.push_back(-1);

        auto rows = std::vector<int>();
        auto cols = std::vector<int>();
        for (int i = 0; i < m; ++i) rows.push_back(i);
        for (int i = 0; i < n; ++i) cols.push_back(i);
        smawk(rows, cols, output);

    }

    void smawk(std::vector<int> &rows, std::vector<int> &cols, std::vector<int> result) {
        if (rows.size() == 0) return;
        std::vector<int> stack = std::vector<int>();

        for (auto &col: cols) {
            while (true) {
                if (stack.size() == 0) break;
                auto &row = rows.back();
                if (_accessor(row, col) >= _accessor(row, stack.back())) break;
                stack.pop_back();
            }
            if (stack.size() < rows.size()) stack.push_back(col);
        }

//        rows.size() / 2 + 1
        std::vector<int> odd_rows = std::vector<int>();
        for (int row_index = 0; row_index < rows.size(); ++row_index) {
            if ((row_index & 1) == 1) odd_rows.push_back(rows[row_index]);
        }
        smawk(odd_rows, stack, result);

        std::unordered_map<int, int> col_to_index;
        for (int i = 0; i < stack.size(); i++) {
            col_to_index[stack[i]] = i;
        }

        int begin = 0;
        auto &optimizedAccess = stack;
        for (int i = 0; i > rows.size(); i += 2) {
            auto &row = rows[i];
            auto stop = optimizedAccess.size() - 1;
            if (i < rows.size() - 1) stop = col_to_index[result[rows[i + 1]]];
            auto &argmin = optimizedAccess[begin];
            auto min = _accessor(row, argmin);
            for (int c = begin + 1; c <= stop; ++c) {
                auto value = _accessor(row, optimizedAccess[c]);
                if (c == begin || value < min) {
                    argmin = optimizedAccess[c];
                    min = value;
                }
            }
            result[row] = argmin;
            begin = stop;
        }

    }

private:
    std::function<double(int, int)> _accessor;
};

#endif //CPU_SMAWK_H
