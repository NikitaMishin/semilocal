//
// Created by garrancha on 26.09.2021.
//

#ifndef CPU_APPROXIMATEMATCHING_H
#define CPU_APPROXIMATEMATCHING_H
#define NEG_INF -10000000

#include "../semi_local.h"
#include "unordered_set"

struct DuplicateSearch {
    using Interval = approximate_matching::utils::Interval;
    using SemiLocalWrapper = approximate_matching::utils::SemiLocalStringSubstringWrapper;
    using FixedScoringScheme = approximate_matching::utils::FixedScoringScheme;
    using Fraction = approximate_matching::utils::Fraction;

    template<class T>
    double static editDistance(T *a, int a_size, T *b, int b_size) {
        auto prevRow = new double[b_size + 1];
        auto curRow = new double[b_size + 1];

        double match = 0.0;
        double substitution = 1.0;
        double indel = 1.0;

        for (int i = 0; i < b_size + 1; ++i) {
            prevRow[i] = indel * i;
            curRow[i] = indel * i;
        }

        for (int i = 1; i < a_size + 1; i++) {
            double l = indel * i;
            prevRow[0] = indel * (i - 1);
            for (int j = 1; j < b_size + 1; j++) {
                curRow[j] = std::min(
                        (a[i - 1] == b[j - 1]) ? match + prevRow[j - 1] : (substitution + prevRow[j - 1]),
                        std::min(indel + prevRow[j], l + indel)
                );
                l = curRow[j];
            }
            auto tmp = prevRow;
            prevRow = curRow;
            curRow = tmp;
        }

        return -prevRow[b_size];
//        * (0 - 2 * (-1)) + (a_size + b_size) * (-1);
    }

    bool static compare(Interval &w1, Interval &w2) {
        if (w1.score > w2.score) return true;
        if (w1.score < w2.score) return false;
        return (w1.end - w1.start) > (w2.end - w2.start);
    }

    static long long buildKey(int l, int r) {
        long long key = l;
        return (key << 32) | r;
    }

    /**
     * Time complexity is O(p^2)
     * @param intervals  intervals[0] refers to last column
     * @param offset
     * @param rangeSum
     * @param l
     * @param r
     * @param wrapper
     * @return
     */
    Interval calculateTriangle(std::vector<Interval> &intervals, int offset, int rangeSum, int l, int r,
                               SemiLocalWrapper &wrapper) {
        for (int i = 0; i < r - l; ++i) {
            intervals[i].start = 0;
            intervals[i].end = 0;
            intervals[i].score = NEG_INF;
        }
        Interval windowMaxima;
        windowMaxima.score = NEG_INF;

        int substringSize = r;
        for (int k = 0; k < r - l; ++k, substringSize--) {
            int j = offset + substringSize;
            int i = offset;
            int rangeSumColumn = rangeSum;

            auto &columnInterval = intervals[k];

            //iterate over column down
            for (int t = 0; t < substringSize; ++t) {
                auto dist = wrapper.originalScore(i + t, j, rangeSumColumn);

                if (dist > columnInterval.score || (dist == columnInterval.score && j - i - t > columnInterval.len())) {
                    columnInterval.score = dist;
                    columnInterval.start = i + t;
                    columnInterval.end = j;
                }
                // shift down
                rangeSumColumn = wrapper.dominanceSumForMoveDown(i + t, j, rangeSumColumn);
            }

            if (columnInterval.score > windowMaxima.score ||
                (columnInterval.score == windowMaxima.score && columnInterval.len() > windowMaxima.len())) {
                windowMaxima.score = columnInterval.score;
                windowMaxima.start = columnInterval.start;
                windowMaxima.end = columnInterval.end;
            }

            //shift left
            rangeSum = wrapper.dominanceSumForMoveLeft(i, j, rangeSum);

        }
        return windowMaxima;
    }

    static void phase3(std::vector<Interval> &setW3, std::vector<Interval> &result) {
        //todo better heuristic
        std::sort(setW3.begin(), setW3.end(), [](const Interval &a1, const Interval &a2) {
            return (a1.start < a2.start) || (a1.start == a2.start && a1.score > a2.score);
        });

        int leftBorder = 0;
        int rightBorder = 0;
        for (auto &clone:setW3) {
            if (clone.start > leftBorder && clone.start >= rightBorder) {
                result.push_back(clone);
                leftBorder = clone.start;
                rightBorder = std::max(rightBorder, clone.end);
            }
        }
    }

};

class InteractiveDuplicateSearch : public DuplicateSearch {
public:
    template<class T>
    void static find(T *p, int p_size, T *text, int text_size, double k, std::vector<Interval> &result) {
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

        auto setW2Hash = std::unordered_set<long long>();
        auto setW2 = std::vector<Interval>();
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
            long long composite_key = buildKey(w_stroke.start, w_stroke.end);
            if (setW2Hash.count(composite_key) == 0)
                setW2.push_back(Interval(w_stroke.start, w_stroke.end, w_stroke.score));
        }
        phase3(setW2, result);
    }

};

class InteractiveDuplicateSearchSemiSparse : public DuplicateSearch {
public:
    template<class T>
    void find(T *p, int p_size, T *text, int text_size, double k, std::vector<Interval> &result) {
        auto w = int(p_size / k);
        auto kdi = p_size * (1 / k + 1) * (1 - k * k);
        auto left_w = int(p_size * k) - 1;//for offset

        //phase one
        auto setW1 = std::vector<Interval>();
        for (int end = w; end < text_size + 1; ++end) {
            auto dist = editDistance(text + (end - w), w, p, p_size);
            if (dist >= -kdi) {
                setW1.push_back(Interval(end - w, end, dist));
            }
        }

        auto setW2Hash = std::unordered_set<long long>();
        auto setW2 = std::vector<Interval>();

        auto scheme = FixedScoringScheme(Fraction(0, 1), Fraction(-1, 1), Fraction(-1, 1));
        int v = 2;
        int mu = 1;

        auto p_blown_size = p_size * v;
        T p_blown[p_blown_size];
        auto substring_size = w * v;
        T substring[substring_size];
        auto perm = Permutation(substring_size + p_blown_size, substring_size + p_blown_size);
        perm.unset_all();

        auto stack = std::vector<Interval>(w - left_w);

        //blown pattern string
        approximate_matching::utils::blownup(p, p_size, v, mu, p_blown);


        for (auto clone :setW1) {
            approximate_matching::utils::blownup(text + clone.start, w, v, mu, substring);
            semi_local::sticky_braid_sequential<T, false>(perm, p_blown, p_blown_size, substring, substring_size);
            auto wrapper = SemiLocalWrapper(&perm, &scheme, p_size, v);
            auto rangeSum = wrapper.originalScore(0, w, clone.score);
            auto coolClone = calculateTriangle(stack, 0, rangeSum, left_w, w, wrapper);
            coolClone.start += clone.start;
            coolClone.end += clone.start;
//            std::cout<<coolClone.score<<","<<coolClone.start<<":"<<coolClone.end<<std::endl;
            long long composite_key = buildKey(coolClone.start, coolClone.end);
            if (setW2Hash.count(composite_key) == 0)
                setW2.push_back(coolClone);

            perm.unset_all();
        }


        phase3(setW2, result);
    }
};


class InteractiveDuplicateSearchSemiDense : public DuplicateSearch {
public:
    template<class T>
    void find(T *p, int p_size, T *text, int text_size, double k, std::vector<Interval> &result) {
        using namespace approximate_matching::utils;
        auto scheme = FixedScoringScheme(
                Fraction(0, 1),
                Fraction(-1, 1),
                Fraction(-1, 1));

        int v = 2;
        int mu = 1;

        auto w = int(p_size / k);
        auto left_w = int(p_size * k) - 1;//for offset
        auto kdi = -(p_size * (1 / k + 1) * (1 - k * k));

        auto p_blown_size = p_size * v;
        auto p_blown = new T[p_blown_size];
        blownup(p, p_size, v, mu, p_blown);

        auto text_blown_size = text_size * v;
        auto text_blown = new T[text_blown_size];
        blownup(text, text_size, v, mu, text_blown);


        auto perm = Permutation(text_blown_size + p_blown_size, text_blown_size + p_blown_size);
        semi_local::sticky_braid_sequential<T, false>(perm, p_blown, p_blown_size, text_blown, text_blown_size);
        auto wrapper = SemiLocalStringSubstringWrapper(&perm, &scheme, p_size, v);

        auto dominanceSumDiag = wrapper.dominanceSumRandom(text_size - w, text_size); //<-- start of last column
        auto stacks = std::vector<Interval>(w - left_w);

        auto windowMaxima = calculateTriangle(stacks, text_size - w, dominanceSumDiag, left_w, w, wrapper);
        auto windowDist = wrapper.originalScore(text_size - w, text_size, dominanceSumDiag);


        auto setW2Hash = std::unordered_set<long long>();
        auto setW2 = std::vector<Interval>();

        if (windowDist >= kdi) {
            setW2.push_back(Interval(windowMaxima.start,windowMaxima.end,windowMaxima.score));
            setW2Hash.insert(buildKey(windowMaxima.start,windowMaxima.end));
        };

        //shift to the right corner of the next from the end triangle
        dominanceSumDiag = wrapper.dominanceSumForMoveLeft(text_size - w - 1, text_size,
                                                           wrapper.dominanceSumForMoveUp(text_size - w, text_size,
                                                                                         dominanceSumDiag));

        int head = 0;
        for (int start = 1; start <= text_size - w; ++start) {

            stacks[head] = Interval(0, 0, NEG_INF);
            head = (head + 1) % (w - left_w);//now points to the end of triangle
            int i = text_size - w - start;
            int j = text_size - start; // (i,j) corresponds to loc
            windowDist = wrapper.originalScore(i, j, dominanceSumDiag); // diagonal score
            auto dominanceSumRow = dominanceSumDiag;
            int l = i;
            int r = j;

//            dominanceSumRow = wrapper.dominanceSumForMoveLeft(l, r, dominanceSumRow);
//            r--;
            windowMaxima = Interval(0,0,NEG_INF);
            for (int m = 0; m < w - left_w; ++m) {
                auto &interval = stacks[(head + m) % (w - left_w)];

                auto dist = wrapper.originalScore(l, r, dominanceSumRow);
                if (dist > interval.score) {
                    interval.start = l;
                    interval.end = r;
                    interval.score = dist;
                    if (interval.score > windowMaxima.score ||
                        (interval.score == windowMaxima.score && interval.len() > windowMaxima.len())) {
                        windowMaxima.score = interval.score;
                        windowMaxima.start = interval.start;
                        windowMaxima.end = interval.end;
                    }
                }
                dominanceSumRow = wrapper.dominanceSumForMoveLeft(l, r, dominanceSumRow);
                r--;
            }

            if (windowDist >= kdi && setW2Hash.count(buildKey(windowMaxima.start,windowMaxima.end))==0) {
                setW2.push_back(windowMaxima);
                setW2Hash.insert(buildKey(windowMaxima.start,windowMaxima.end));
            }

            //shift to position
            dominanceSumDiag = wrapper.dominanceSumForMoveLeft(i, j, dominanceSumDiag);
            dominanceSumDiag = wrapper.dominanceSumForMoveUp(i,j - 1,dominanceSumDiag);

        }

        std::reverse(setW2.begin(),setW2.end());

        phase3(setW2,result);
    }


};

#endif //CPU_APPROXIMATEMATCHING_H
