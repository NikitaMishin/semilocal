#ifndef CPU_LIBRARY_H
#define CPU_LIBRARY_H

#include <vector>


#include <iostream>
#include <bitset>
#include <cstring>

/**
 *
 * @tparam Input
 * @param a
 * @param b
 * @return
 */
template<class Input>
int naive_prefix_lcs(std::vector<Input> a, std::vector<Input> b) {
    int arr[a.size() + 1][b.size() + 1];
    auto m = a.size() + 1;
    auto n = b.size() + 1;
    for (auto j = 0; j < n; j++) {
        arr[0][j] = 0;
    }

    for (auto i = 0; i < m; i++) {
        arr[i][0] = 0;
    }
    for (int i = 1; i < m; ++i) {
        for (int j = 1; j < n; ++j) {
            arr[i][j] = std::max(std::max(arr[i - 1][j], arr[i][j - 1]),
                                 (a[i - 1] == b[j - 1]) ? arr[i - 1][j - 1] + 1 : arr[i - 1][j - 1]);
        }
    }

//    for (int i = 0; i < m; ++i) {
//        for (int j = 0; j < n; ++j) {
//            std::cout<<arr[i][j] <<" ";
//        }
//        std::cout<<std::endl;
//    }

    return arr[m - 1][n - 1];
}

/**
 *
 * @tparam Input
 * @param a
 * @param b
 * @return
 */
template<class Input>
int prefix_lcs_sequential(std::vector<Input> a, std::vector<Input> b) {

    std::vector<Input> input_a;
    std::vector<Input> input_b;
    int m, n;

    if (a.size() > b.size()) {
        m = a.size() + 1;
        n = b.size() + 1;
        input_a = a;
        input_b = b;
    } else {
        n = a.size() + 1;
        m = b.size() + 1;
        input_b = a;
        input_a = b;
    }

    auto prev_row = new int[n];
    auto cur_row = new int[n];
    for (int i = 0; i < n; ++i) {
        cur_row[i] = 0;
        prev_row[i] = 0;
    }

    for (int i = 1; i < m; ++i) {
        auto l = 0;
        for (int j = 1; j < n; ++j) {
            cur_row[j] = std::max(
                    std::max(prev_row[j], l),
                    (input_a[i - 1] == input_b[j - 1]) ? prev_row[j - 1] + 1 : prev_row[j - 1]
            );
            l = cur_row[j];

        }
        std::swap(prev_row, cur_row);
    }

    return prev_row[n - 1];

}


template<class Input>
int prefix_lcs_sequential_skewed(std::vector<Input> a, std::vector<Input> b) {
//    check special case 2x2

    if (a.size() == 1 && b.size() == 1) {
        return a[0] == b[0] ? 1 : 0;
    }

    auto diagonal_size = 1 + std::min(a.size(), b.size());
    auto a1 = new int[diagonal_size];
    auto a2 = new int[diagonal_size];
    auto a3 = new int[diagonal_size];

    auto pos_i = 0;
    auto pos_j = 0;
    auto start_i = pos_i;
    auto start_j = pos_j;
    auto min = std::min(a.size(), b.size());
    auto num_diag = a.size() + b.size() + 1;
    auto total_same_length_diag = num_diag - (min + 1) - min;

//    init step
    for (int k = 0; k < diagonal_size; ++k) {
        a3[k] = 0;
        a2[k] = 0;
        a1[k] = 0;
    }


    start_i--;
    // fill upper square
    for (int k = 2; k <= diagonal_size; ++k, start_i++) {
        pos_i = start_i;
        pos_j = 0;
        a3[0] = 0;
        a3[k - 1] = 0;
        for (int i = 1; i < k - 1; ++i) {
            a3[i] = std::max(
                    std::max(a2[i], a2[i - 1]),
                    (a[pos_i] == b[pos_j]) ? 1 + a1[i - 1] : a1[i - 1]
            );
            pos_i--;
            pos_j++;
        }
        std::swap(a1, a2);
        std::swap(a3, a2);
    }
//    if (a.size() <= b.size()) {
//        start_i--;
//    }

    // phase 2:: fill
    if (a.size() >= b.size()) {
        //        same pattern
//        todo length
        for (int k = 0; k < total_same_length_diag; ++k, start_i++) {
            pos_i = start_i;
            pos_j = 0;
            a3[0] = 0;

            for (int i = 1; i < diagonal_size; ++i) {
                a3[i] = std::max(
                        std::max(a2[i], a2[i - 1]),
                        (a[pos_i] == b[pos_j]) ? 1 + a1[i - 1] : a1[i - 1]
                );
                pos_i--;
                pos_j++;
            }
            std::swap(a1, a2);
            std::swap(a3, a2);
        }

    }

    //      special case when:
    //      a==b => |a1| = c-1 , |a2| = c, |a3|= c-1  or
    //      a>b  => |a1| = c, |a2| = c, |a3| = c-1
    //      a<b ->  |a1| = c - 1, |a2| = c, |a3| = c

    pos_i = start_i;
    pos_j = 0;

    if (a.size() < b.size()) {
//        todo
        a3[diagonal_size - 1] = 0;
//        total_same_length_diag--; //since we fill it now
    }

    for (int i = 0; i < diagonal_size - 1; ++i) {
        a3[i] = std::max(
                std::max(a2[i], a2[i + 1]),
                (a[pos_i] == b[pos_j]) ? 1 + a1[i] : a1[i]
        );
        pos_i--;
        pos_j++;
    }
    start_j++;
    std::swap(a1, a2);
    std::swap(a3, a2);

    if (a.size() < b.size()) {
//        since special case then -1
        for (int k = 0; k < total_same_length_diag; ++k, start_j++) {
            pos_i = start_i;
            pos_j = start_j;

            a3[diagonal_size - 1] = 0;
            for (int i = 0; i < diagonal_size - 1; ++i) {
                a3[i] = std::max(
                        std::max(a2[i], a2[i + 1]),
                        (a[pos_i] == b[pos_j]) ? 1 + a1[i + 1] : a1[i + 1]
                );
                pos_i--;
                pos_j++;
            }
            std::swap(a1, a2);
            std::swap(a3, a2);
        }
    }

    if (a.size() >= b.size()) diagonal_size -= 1;


//    phase 3
//    pattern a3[i] = max(a1[i+1],a2[i],a2[i-1])
    for (int size = diagonal_size - 1; size > 1; size--, start_j++) {
        pos_i = start_i;
        pos_j = start_j;

        for (int i = 0; i < size; ++i) {
            a3[i] = std::max(
                    std::max(a2[i], a2[i + 1]),
                    (a[pos_i] == b[pos_j]) ? 1 + a1[i + 1] : a1[i + 1]
            );
            pos_i--;
            pos_j++;
        }

        std::swap(a1, a2);
        std::swap(a3, a2);
    }


    //  need to calculate last one cell
    return std::max(std::max(a2[0], a2[1]), (a[a.size() - 1]) == b[b.size() - 1] ? 1 + a1[1] : a1[1]);
}


std::vector<int> gen_seq(int size, int alphabet_size = 5) {

    auto v = std::vector<int>();
    for (int i = 0; i < size; ++i) {
        v.push_back(rand() % alphabet_size);
    }
    return v;
}

/**
 * Compute via antidiagonal pattern
 * @tparam Input
 * @param a
 * @param b
 * @return
 */
template<class Input>
int prefix_lcs_parallel_mpi(std::vector<Input> a, std::vector<Input> b, int threads_num = 2) {
    auto max_diagonal_size = 1 + std::min(a.size(), b.size());
    auto diagonal_prev = new int[max_diagonal_size];
    auto diagonal_cur = new int[max_diagonal_size];
    auto num_diagonals = a.size() + b.size() + 1;
//    fill with zeros with threads

//

}


template<class Input>
int prefix_lcs_via_braid_sequential(std::vector<Input> a, std::vector<Input> b) {

    auto m = a.size();
    auto n = b.size();

    auto dis_braid = 0;
    auto size = m + n;
    auto strand_map  = new bool[size];

    for (int k = 0; k < m; ++k) {
        strand_map[k] = true;
    }

    for (int l = m; l < m + n; ++l) {
        strand_map[l] = false;
    }


    for (int i = 0, left_edge = m - 1; i < m; ++i, left_edge--) {
        for (int j = 0, top_edge = m; j < n; ++j, top_edge++) {
            bool left_strand =  strand_map[left_edge];
            bool right_strand = strand_map[top_edge];
            if ( a[i] == b[j] || (!left_strand && right_strand)) {
                strand_map[left_edge] = right_strand;
                strand_map[top_edge] = left_strand;
            }

        }
    }


    for (int i1 = 0; i1 < m; ++i1) {
        dis_braid+=strand_map[i1];
    }


    return   m - dis_braid;
}


template<class Input>
int sticky_braid_sequential(std::vector<Input> a, std::vector<Input> b) {
    auto m = a.size();
    auto n = b.size();

    auto size = m + n;
    auto strand_map = new int[size];
    auto reduced_sticky_braid = new int[size];
    for (int i = 0; i < size; ++i) {
        strand_map[i] = i;
    }
    auto ss =0;

    for (int i = 0; i < m; ++i) {
        for (int j = 0; j < n; ++j) {
            auto left_edge = m - 1 - i;
            auto top_edge = m + j;
            auto left_strand = strand_map[left_edge];
            auto right_strand = strand_map[top_edge];

            if (a[i] == b[j] || (a[i] != b[j] && left_strand > right_strand)) {
                strand_map[left_edge] = right_strand;
                strand_map[top_edge] = left_strand;
            }

            if (j == n - 1) {
                auto strand_end = left_edge + n;
                auto strand_start = strand_map[left_edge];
                reduced_sticky_braid[strand_start] = strand_end;
                if (strand_start < m) ss++;
            }

            if (i == m - 1) {
                auto strand_end = top_edge - m;
                auto strand_start = strand_map[top_edge];
                reduced_sticky_braid[strand_start] = strand_end;
            }
        }
    }
    return m-ss;
}


#endif //CPU_LIBRARY_H
