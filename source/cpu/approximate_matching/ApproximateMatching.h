//
// Created by garrancha on 26.09.2021.
//

#ifndef CPU_APPROXIMATEMATCHING_H
#define CPU_APPROXIMATEMATCHING_H
#define NEG_INF -10000000

#include "../semi_local.h"

class InteractiveDuplicateSearch {
public:
    template<class T>
    void static find(T *p, int p_size, T *text, int text_size, double k,
                     std::vector<approximate_matching::utils::Interval> &result) {
        using namespace approximate_matching::utils;
        auto w = int(p_size / k);
        auto kdi = p_size * (1 / k + 1) * (1 - k * k);


        //phase one
        auto setW1 = std::vector<Interval>();
        for (int end = w; end < text_size + 1; ++end) {
            auto dist = editDistance(text + (end - w), w, p, p_size);
            if (dist >= -kdi) {
                setW1.push_back(Interval(end - w, end, dist));
            }
        }

        auto setW2 = std::map<long long, Interval>();
        //phase 2
        for (auto &clone :setW1) {
            auto w_stroke = clone;
            // iterated over all sizes
            for (int l = int(p_size * k); l <= w; l++) {
                for (int end = l + clone.start; end <= clone.end; end++) {
                    auto w2 = Interval(end - l, end, editDistance(p, p_size, text + (end - l), l));
                    if (compare(w2, w_stroke)) {
                        w_stroke = w2;
                    }
                }
            }
            long long composite_key = w_stroke.start;
            composite_key = (composite_key << 32) + w_stroke.end;
            if (setW2.count(composite_key) == 0)
                setW2.insert({composite_key, Interval(w_stroke.start, w_stroke.end, w_stroke.score)});
        }


//        //last comparator refers big one or a small one
        std::vector<Interval> setW3;
        for (auto &clone: setW2) setW3.push_back(clone.second);

        //todo better heuristic
        std::sort(setW3.begin(), setW3.end(), [](const Interval &a1, const Interval &a2) {
            return (a1.start < a2.start) || (a1.start == a2.start && a1.score > a2.score);
        });

        int leftBorder = 0;
        int rightBorder = 0;
        for (auto &clone:setW3) {
//            std::cout<<clone.start<<":"<<clone.end<<"score="<<clone.score<<std::endl;
            if (clone.start > leftBorder && clone.start >= rightBorder) {
//                std::cout << clone.start << ":" << clone.end << std::endl;

                result.push_back(clone);
                leftBorder = clone.start;
                rightBorder = std::max(rightBorder, clone.end);
            }
        }


    }

//private:
    template<class T>
    double static editDistance(T *a, int a_size, T *b, int b_size) {
        auto prevRow = new double[b_size + 1];
        auto curRow = new double[b_size + 1];
        for (int i = 0; i < b_size + 1; ++i) {
            prevRow[i] = 0.0;
            curRow[i] = 0.0;
        }
        double match = 1.0;
        double mismatch = 0.0;
        double gap = 0.0;

        for (int i = 1; i < a_size + 1; i++) {
            double l = 0.0;
            for (int j = 1; j < b_size + 1; j++) {
                curRow[j] = std::max(
                        (a[i - 1] == b[j - 1]) ? match + prevRow[j - 1] : mismatch + prevRow[j - 1],
                        std::max(gap + prevRow[j], l + gap)
                );
                l = curRow[j];
            }
            auto tmp = prevRow;
            prevRow = curRow;
            curRow = tmp;
        }

        return prevRow[b_size] * (0 - 2 * (-1)) + (a_size + b_size) * (-1);
    }

    bool static compare(approximate_matching::utils::Interval &w1, approximate_matching::utils::Interval &w2) {
        if (w1.score > w2.score) return true;
        if (w1.score < w2.score) return false;
        return (w1.end - w1.start) > (w2.end - w2.start);
    }
};


class InteractiveDuplicateSearchSemi {
public:
    template<class T>
    void static find(T *p, int p_size, T *text, int text_size, double k,
                     std::vector<approximate_matching::utils::Interval> &result) {
        using namespace approximate_matching::utils;
        auto scheme = FixedScoringScheme(Fraction(0, 1), Fraction(-1, 1), Fraction(-1, 1));

        int v = 2;
        int mu = 1;

        auto p_blown_size = p_size * v;
        auto p_blown = new T[p_blown_size];

        blownup(p, p_size, v, mu, p_blown);

        auto text_blown_size = text_size * v;
        auto text_blown = new T[text_blown_size];
        blownup(text, text_size, v, mu, text_blown);


        auto perm = Permutation(text_blown_size + p_blown_size, text_blown_size + p_blown_size);


        auto w = int(p_size / k);
        auto left_w = int(p_size * k) - 1;//for offset
        auto kdi = p_size * (1 / k + 1) * (1 - k * k);


        semi_local::sticky_braid_sequential<T, false>(perm, p_blown, p_blown_size, text_blown, text_blown_size);

        auto wrapper = SemiLocalStringSubstringWrapper(&perm,&scheme,p_size,v);

        auto dominanceSumDiag = wrapper.dominanceSumRandom(text_size - w, text_size); //<-- start of last column

        auto stacks = std::vector<Interval>(w - left_w);

        auto windowMaxima = calculateTriangle(stacks, text_size - w, dominanceSumDiag, left_w, w, wrapper);

        auto windowDist = wrapper.originalScore(text_size - w, text_size, dominanceSumDiag);

        std::map<long long,Interval> setW2;


        if(windowDist >= -kdi)setW2.insert({buildKey(windowMaxima.start,windowMaxima.end),windowMaxima});
        dominanceSumDiag = wrapper.dominanceSumForMoveLeft(text_size - w - 1,text_size,
                                                           wrapper.dominanceSumForMoveUp(text_size-w,text_size,dominanceSumDiag));
        int head = 1;

        for (int start = 1; start <= text_size - w; ++start) {
            head = head % (w - left_w);
            stacks[head] = Interval(0,0,NEG_INF);
            int i = text_size - w - start;
            int j = text_size - start;
            windowDist = wrapper.originalScore(i, j, dominanceSumDiag);

            auto dominanceSumRow = dominanceSumDiag;
            int l = i;
            int r = j;

            dominanceSumRow = wrapper.dominanceSumForMoveLeft(l, r, dominanceSumRow);
            r--;

            for (int m = 1; m < w - left_w; ++m) {
                auto &interval = stacks[(head + m) % (w - left_w)];

                auto dist = wrapper.originalScore(l, r, dominanceSumRow);
                if (dist > interval.score) {
                    interval.start = l;
                    interval.end = r;
                    interval.score = dist;
                    if (interval.score > localMaxima.score ||
                        (interval.score == localMaxima.score && interval.len() > localMaxima.len())) {
                        localMaxima.score = interval.score;
                        localMaxima.start = interval.start;
                        localMaxima.end = interval.end;
                    }
                }
                dominanceSumRow = wrapper.dominanceSumForMoveLeft(l, r, dominanceSumRow);
                r--;
            }

            if (windowDist >= -kdi) {
                //todo add to set localMaxima
            }


            dominanceSumDiag = dominanceSumForMoveLeft(i, j, dominanceSumDiag, v);
            j--;
            dominanceSumDiag = dominanceSumForMoveUp(i, j, dominanceSumDiag, v);
            head = (head + 1) % w;
//            todo clear
//            stakc
        }

        // calculate for a last square
        auto dominanceSumStart = dominanceSumRandom(20, 27);
        std::cout << "dominance sum:" << dominanceSumStart << std::endl;
        auto value = canonicalDecompositionWithKnown(20, 27, v, p_size, dominanceSumStart);
        std::cout << "value:" << value << std::endl;
        std::cout << "reg:" << originalScore(20, 27, dominanceSumStart, v, p_size) << std::endl;



//        std::cout<<canonicalDecomposition(text_size - p_size + p_size, text_size, v, p_size);
//        size from 1 to 3
//   0  1 2 3 4    [1,2]                              (0,1),(0,2),(0,3)
//           +     +      +                                 col1  col2 col3
// 0 (0,0) (0,1) (0,2), (0,3), (0,4)                        upd   upd   del
// 1     - (1,1) (1,2), (1,3), (1,4) +                     (1,2),(1,3),(1,4) <---here
// 2 - -         (2,2)  (2,3), (2,4) +                           (2,3),(2,4)  <--- window
// 3 - - -              (3,3)  (3,4) + <                               (3,4)
// 4 - - - -                   (4,4)
//
// 1,4 вниз 3,
// 0,3 вниз 3
//  обработали (3,4)
//  0,3 all

//        //last comparator refers big one or a small one
//        std::vector<Interval> setW3;
//        for (auto &clone: setW2) setW3.push_back(clone.second);
//
//        //todo better heuristic
//        std::sort(setW3.begin(), setW3.end(), [](const Interval &a1, const Interval &a2) {
//            return (a1.start < a2.start) || (a1.start == a2.start && a1.score > a2.score);
//        });
//
//        int leftBorder = 0;
//        int rightBorder = 0;
//        for (auto &clone:setW3) {
////            std::cout<<clone.start<<":"<<clone.end<<"score="<<clone.score<<std::endl;
//            if (clone.start > leftBorder && clone.start >= rightBorder) {
////                std::cout << clone.start << ":" << clone.end << std::endl;
//
//                result.push_back(clone);
//                leftBorder = clone.start;
//                rightBorder = std::max(rightBorder, clone.end);
//            }
//        }


    }

    static long long buildKey(int l,int r){
        long long key = l;
        return  (key<<32) | r;
    }
    approximate_matching::utils::Interval
    static calculateTriangle(std::vector<approximate_matching::utils::Interval>&intervals,int offset,int rangeSum,int l,int r, approximate_matching::utils::SemiLocalStringSubstringWrapper& wrapper){
        for (int i = 0; i < r - l ; ++i) {
            intervals[i].start = 0;
            intervals[i].end = 0;
            intervals[i].score = NEG_INF;
        }
        approximate_matching::utils::Interval windowMaxima;
        windowMaxima.score = NEG_INF;

        int substringSize = r;
        for (int k = 0; k <  r - l; ++k, substringSize--) {
            int j = offset + substringSize;
            int i = offset;
            int rangeSumColumn = rangeSum;
            //shift left
            rangeSum = wrapper.dominanceSumForMoveLeft(i, j, rangeSum);

            auto& columnInterval = intervals[k];

            for (int t = 0; t < substringSize; ++t) {
                auto dist = wrapper.originalScore(i + t, j, rangeSumColumn);

                if(dist > columnInterval.score ||(dist==columnInterval.score && j - i - t  > columnInterval.len())) {
                    columnInterval.score = dist;
                    columnInterval.start = i + t;
                    columnInterval.end = j;
                }
                // shift down
                wrapper.dominanceSumForMoveDown(i + t, j, rangeSumColumn);
            }

            if (columnInterval.score > windowMaxima.score || (columnInterval.score==windowMaxima.score&& columnInterval.len() > windowMaxima.len())){
                windowMaxima.score = columnInterval.score;
                windowMaxima.start = columnInterval.start;
                windowMaxima.end = columnInterval.end;
            }
        }
        return windowMaxima;
    }


};

#endif //CPU_APPROXIMATEMATCHING_H
