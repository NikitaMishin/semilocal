//
// Created by garrancha on 26.09.2021.
//

#ifndef CPU_APPROXIMATEMATCHING_H
#define CPU_APPROXIMATEMATCHING_H
#define NEG_INF (-10000000)

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


        auto value = -prevRow[b_size];

        delete[] prevRow;
        delete[] curRow;
        return value;
    }


    template<class T>
    Interval findClosestLongestSubstring(T *p, int p_size, T *substring, int min_size, int max_size, int offset,
                                         std::unordered_map<long long, double> &precalc, double kdi) {
        Interval interval = Interval(0, 0, NEG_INF);
/*
        //iterative over row backward starting from largest
        for (int length = max_size; length >= min_size ; length--) {
            int j = length;
            //iterate  over column from top to bottom
            for (int i = 0; i <= length - min_size; ++i) {
                double edit;
                auto key = buildKey(i+offset,offset + j);
                if(precalc.count(key)){
                    edit = precalc[key];
                } else {
                    edit = editDistance(p,p_size, substring + i, j - i);
                    precalc[key] = edit;
                }
                auto potential = Interval(i + offset, offset + j, edit);
                if(compare(potential,interval)){
                    interval = Interval(potential.start,potential.end,potential.score);
                }
            }
        }*/

///*
        for (int length = min_size; length <= max_size; length++) {
            for (int end = length; end < max_size + 1;) {
                double edit;
                int i = end - length;
                int j = end;
                auto key = buildKey(i + offset, offset + j);
                if (precalc.count(key)) {
                    edit = precalc[key];
                } else {
                    edit = editDistance(p, p_size, substring + i, j -  i);
                    precalc[key] = edit;
                }
                auto potential = Interval(i + offset, offset + j, edit);
                if (compareUPD(potential, interval)) {
                    interval = Interval(potential.start, potential.end, potential.score);
                }
                if (edit < interval.score - 1) {
                    end += std::max(1, int((abs(edit) - abs(interval.score)) / 2));
                } else {
                    ++end;
                }
            }

        }
//        */


        return interval;
    }

    bool static compare(Interval &w1, Interval &w2) {
        if (w1.score > w2.score) return true;
        if (w1.score < w2.score) return false;
        return (w1.len()) > (w2.len());
    }

    bool static compareUPD(Interval &w1, Interval &w2) {
        if (w1.score > w2.score) return true;
        if (w1.score < w2.score) return false;
        return (w1.len()) >= (w2.len());
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
        //iterate over row left
        for (int k = 0; k < r - l; ++k, substringSize--) {
            int i = offset;
            int j = offset + substringSize;
            int rangeSumColumn = rangeSum;

            auto &columnInterval = intervals[k];

            //iterate over column down
            for (int t = 0; t < substringSize - l; ++t) {
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
    void find(T *p, int p_size, T *text, int text_size, double k, std::vector<Interval> &result) {
        auto w = int(p_size / k);
        auto kdi = -(p_size * (1 / k + 1) * (1 - k * k));
        auto left_w = int(p_size * k);

        auto precalc = std::unordered_map<long long, double>();


        //phase one
        auto setW1 = std::vector<Interval>();
        for (int end = w; end < text_size + 1;) {
            auto dist = editDistance(text + (end - w), w, p, p_size);
            precalc[buildKey(end - w, end)] = dist;
            if (dist >= kdi) setW1.push_back(Interval(end - w, end, dist));

            if (dist < kdi - 1) {
                end += std::max(1, int((abs(dist) - abs(kdi)) / 2));
            } else {
                ++end;
            }
        }

        auto setW2Hash = std::unordered_set<long long>();
        auto setW2 = std::vector<Interval>();
        //phase 2
        int position = 0;
        for (auto &clone :setW1) {
            auto w_stroke = findClosestLongestSubstring(p, p_size, text + clone.start, left_w, w, clone.start, precalc,
                                                        kdi);
            long long composite_key = buildKey(w_stroke.start, w_stroke.end);
            if (setW2Hash.count(composite_key) == 0) {
                setW2Hash.insert(composite_key);
                setW2.push_back(Interval(w_stroke.start, w_stroke.end, w_stroke.score));
            }

            position++;
        }

        phase3(setW2, result);
    }

};

class InteractiveDuplicateSearchSemiSparse : public DuplicateSearch {
public:
    template<class T>
    void find(T *p, int p_size, T *text, int text_size, double k, std::vector<Interval> &result) {
        auto w = int(p_size / k);
        auto kdi = -(p_size * (1 / k + 1) * (1 - k * k));
        auto left_w = int(p_size * k) -1;//for offset

        //phase one
        auto setW1 = std::vector<Interval>();
        for (int end = w; end < text_size + 1;) {
            auto dist = editDistance(text + (end - w), w, p, p_size);
            if (dist >= kdi) setW1.push_back(Interval(end - w, end, dist));

            if (dist < kdi - 1) {
                end += std::max(1, int((abs(dist) - abs(kdi)) / 2));
            } else {
                ++end;
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


        int position = 0;
        for (auto &clone :setW1) {
            approximate_matching::utils::blownup(text + clone.start, w, v, mu, substring);
            semi_local::sticky_braid_sequential<T, false>(perm, p_blown, p_blown_size, substring, substring_size);
            auto wrapper = SemiLocalWrapper(&perm, &scheme, p_size, v);
            auto rangeSum = wrapper.originalScore(0, w, clone.score);
            auto coolClone = calculateTriangle(stack, 0, rangeSum, left_w, w, wrapper);
            coolClone.start += clone.start;
            coolClone.end += clone.start;
            long long composite_key = buildKey(coolClone.start, coolClone.end);
            if (setW2Hash.count(composite_key) == 0) {
                setW2Hash.insert(composite_key);
                setW2.push_back(coolClone);
            }


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
            setW2.push_back(Interval(windowMaxima.start, windowMaxima.end, windowMaxima.score));
            setW2Hash.insert(buildKey(windowMaxima.start, windowMaxima.end));
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


            windowMaxima = Interval(0, 0, NEG_INF);
            for (int m = 0; m < w - left_w; ++m, r--) {
                auto &interval = stacks[(head + m) % (w - left_w)];
                auto dist = wrapper.originalScore(l, r, dominanceSumRow);
                auto potential = Interval(l, r, dist);
                if (compare(potential, interval)) {
                    interval.start = potential.start;
                    interval.end = potential.end;
                    interval.score = potential.score;
                }
                if (compare(interval, windowMaxima)) {
                    windowMaxima.score = interval.score;
                    windowMaxima.start = interval.start;
                    windowMaxima.end = interval.end;
                }
                dominanceSumRow = wrapper.dominanceSumForMoveLeft(l, r, dominanceSumRow);
            }

            if (windowDist >= kdi && setW2Hash.count(buildKey(windowMaxima.start, windowMaxima.end)) == 0) {
                setW2.push_back(windowMaxima);
                setW2Hash.insert(buildKey(windowMaxima.start, windowMaxima.end));
            }

            //shift to position
            dominanceSumDiag = wrapper.dominanceSumForMoveLeft(i, j, dominanceSumDiag);
            dominanceSumDiag = wrapper.dominanceSumForMoveUp(i, j - 1, dominanceSumDiag);

        }

        std::reverse(setW2.begin(), setW2.end());


        phase3(setW2, result);
        delete[] p_blown;
        delete[] text_blown;
    }


    template<class T>
    void print(T *string, int start, int end) {
        std::cout << "\"";
        for (int i = start; i < end; ++i) {
            std::cout << char(string[i]);
        }
        std::cout << "\"" << std::endl;
    }


};

#endif //CPU_APPROXIMATEMATCHING_H
