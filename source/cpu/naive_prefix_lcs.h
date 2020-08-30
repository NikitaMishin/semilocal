//
// Created by nikita on 30.07.2020.
//

#ifndef CPU_PREFIX_LCS_H
#define CPU_NAIVE_PREFIX_LCS_H

#include <vector>
#include <cmath>
#include "omp.h"

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
        a3[diagonal_size - 1] = 0;
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




//
///**
// * Compute via antidiagonal pattern using open mpi
// * @tparam Input
// * @param a
// * @param b
// * @return
// */
//template<class Input>
//int prefix_lcs_parallel_mpi(std::vector<Input> a, std::vector<Input> b, int threads_num = 2) {
//
//    if (a.size() == 1 && b.size() == 1) {
//        return a[0] == b[0] ? 1 : 0;
//    }
//
//    auto mn = 1 + std::min(a.size(), b.size());
//    auto a1 = new int[mn];
//    auto a2 = new int[mn];
//    auto a3 = new int[mn];
//
////
//#pragma omp parallel num_threads(threads_num) default(none) firstprivate(threads_num) shared(a1, a2, a3, , a, b,std::cout)
//{
//
//
//    auto result = 0;
//    #pragma omp parallel num_threads(threads_num) default(none) firstprivate(threads_num) shared(a1, a2, a3, result, a, b,std::cout)
//    {
//        int diagonal_size = 1 + std::min(a.size(), b.size());
//        auto pos_i = 0;
//        auto pos_j = 0;
//
//        auto start_i = -1;
//        auto start_j = pos_j;
//
//        auto min = std::min(a.size(), b.size());
//        auto num_diag = a.size() + b.size() + 1;
//        auto total_same_length_diag = num_diag - (min + 1) - min;
//
//
//        // create private pointers to int
//        int *a1_ptr = a1;
//        int *a2_ptr = a2;
//        int *a3_ptr = a3;
//
//        auto thread_id = omp_get_thread_num();
////        std::cout<<result;
//
//        // responsible for interevals that current thread proccess at specific diagonal
//        // if out of bound then thread do nothing and just swap pointers a1,a2,a3
//        auto max_diag_elem_per_thread = int(std::ceil((diagonal_size * 1.0) / threads_num));
//        auto upper_bound = std::min(max_diag_elem_per_thread * (thread_id + 1), diagonal_size);
//        auto lower_bound = thread_id * max_diag_elem_per_thread;
//
//        //    init step, fill with zeros
//        for (int k = lower_bound; k < upper_bound; ++k) {
//            a3[k] = 0;
//            a2[k] = 0;
//            a1[k] = 0;
//        }
//
//
//        #pragma omp barrier
//
//        // fill upper square
//        for (int k = 2; k <= diagonal_size; ++k, start_i++) {
//            auto diag_elems_per_thread = int(std::ceil(((k - 2) * 1.0) / threads_num));
//
//            pos_i = start_i - thread_id * diag_elems_per_thread;
//            pos_j = 0 + thread_id * diag_elems_per_thread; //different to each thread
//            a3_ptr[0] = 0;
//            a3_ptr[k - 1] = 0;
//            lower_bound = 1 + (diag_elems_per_thread * thread_id);
//            upper_bound = std::min(1 + (thread_id + 1) * diag_elems_per_thread, k - 2);
//
//            for (int i = lower_bound; i < upper_bound; ++i, pos_i--, pos_j++) {
//                a3[i] = std::max(std::max(a2_ptr[i], a2_ptr[i - 1]),
//                                 (a[pos_i] == b[pos_j]) ? 1 + a1_ptr[i - 1] : a1_ptr[i - 1]);
//            }
//            //each thread swaps pointers
//            std::swap(a1_ptr, a2_ptr);
//            std::swap(a3_ptr, a2_ptr);
//            #pragma omp barrier
//        }
//
//
//
//        // phase 2:: fill
//        if (a.size() >= b.size()) {
//
//            for (int k = 0; k < total_same_length_diag; ++k, start_i++) {
//                auto diag_elems_per_thread = int(std::ceil(((diagonal_size - 1) * 1.0) / threads_num));
//                pos_i = start_i - thread_id * diag_elems_per_thread;
//                pos_j = 0 + thread_id * diag_elems_per_thread;
//                a3_ptr[0] = 0;
//                lower_bound = 1 + (thread_id) * diag_elems_per_thread;
//                upper_bound = std::min(1 + (thread_id + 1) * diag_elems_per_thread, diagonal_size - 1);
//
//                for (int i = lower_bound; i < upper_bound; ++i, pos_i--, pos_j++) {
//                    a3_ptr[i] = std::max(std::max(a2_ptr[i], a2_ptr[i - 1]),
//                                         (a[pos_i] == b[pos_j]) ? 1 + a1_ptr[i - 1] : a1_ptr[i - 1]);
//                }
//
//                std::swap(a1_ptr, a2_ptr);
//                std::swap(a3_ptr, a2_ptr);
//#pragma omp barrier
//            }
//        }
//
//
//        if (a.size() < b.size()) {
//            a3_ptr[diagonal_size - 1] = 0;
//        }
//
////        no need to sync
//        auto per_thread = int(std::ceil(((diagonal_size - 1) * 1.0) / threads_num));
//        lower_bound = thread_id * per_thread;
//        upper_bound = std::min(diagonal_size - 1, (thread_id + 1) * per_thread);
//        pos_i = start_i - thread_id * per_thread;
//        pos_j = per_thread * thread_id;
//
//        for (int i = lower_bound; i < upper_bound; ++i, pos_j++, pos_i--) {
//            a3_ptr[i] = std::max(std::max(a2_ptr[i], a2_ptr[i + 1]),
//                                 (a[pos_i] == b[pos_j]) ? 1 + a1_ptr[i] : a1_ptr[i]);
//        }
//
//        start_j++;
//        std::swap(a1_ptr, a2_ptr);
//        std::swap(a3_ptr, a2_ptr);
//#pragma omp barrier
//
//
//        if (a.size() < b.size()) {
//            per_thread = int(std::ceil(((diagonal_size - 1) * 1.0) / threads_num));
//            lower_bound = thread_id * per_thread;
//            upper_bound = std::min(diagonal_size - 1, (thread_id + 1) * per_thread);
//
//
//            for (int k = 0; k < total_same_length_diag; ++k, start_j++) {
//                pos_i = start_i - thread_id * per_thread;
//                pos_j = start_j + per_thread * thread_id;
////todo
//                a3_ptr[diagonal_size - 1] = 0;
//
//                for (int i = lower_bound; i < upper_bound; ++i, pos_i--, pos_j++) {
//                    a3_ptr[i] = std::max(std::max(a2_ptr[i], a2_ptr[i + 1]),
//                                         (a[pos_i] == b[pos_j]) ? 1 + a1_ptr[i + 1] : a1_ptr[i + 1]);
//                }
//                std::swap(a1_ptr, a2_ptr);
//                std::swap(a3_ptr, a2_ptr);
//        #pragma  omp barrier
//            }
//        }
//
//        if (a.size() >= b.size()) diagonal_size -= 1;
//
//
//        for (int size = diagonal_size - 1; size > 1; size--, start_j++) {
//            per_thread = int(std::ceil(((size) * 1.0) / threads_num));
//
//            pos_i = start_i - per_thread * thread_id;
//            pos_j = start_j + per_thread * thread_id - 1 ;
//            lower_bound = thread_id * per_thread;
//            upper_bound = std::min(size, (thread_id + 1) * per_thread);
//            for (int i = lower_bound; i < upper_bound; ++i) {
////                if (pos_j >=b.size()) std::cout<<"J"<<pos_j<<","<<size<<std::endl;
//                a3[i] = std::max(
//                        std::max(a2[i], a2[i + 1]),
//                        (a[pos_i] == b[pos_j]) ? 1 + a1[i + 1] : a1[i + 1]
//                );
//                pos_i--;
//                pos_j++;
//            }
//
//            std::swap(a1_ptr, a2_ptr);
//            std::swap(a3_ptr, a2_ptr);
//            #pragma omp barrier
//        }
//
//        //master thread
//        #pragma omp master
//        {
//            result = std::max(std::max(a2_ptr[0], a2_ptr[1]),
//                              (a[a.size() - 1]) == b[b.size() - 1] ? 1 + a1_ptr[1] : a1_ptr[1]);
//        }
////#pragma  omp barrier
//
//    }
//    //  need to calculate last one cell
//    return result;
//
//}
//

#endif //CPU_PREFIX_LCS_H
