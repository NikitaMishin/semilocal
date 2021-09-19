//
// Created by garrancha on 13.09.2021.
//

#ifndef CPU_WINDOWAMATCH_H
#define CPU_WINDOWAMATCH_H
#include "utils.h"
#include "../approximate_matching/CompleteAMatch.h"
#include "../unit_monge_mult/matrices.h"
#include "../unit_monge_mult/dominance_sum_queries.h"
template<class Integer>
class IThresholdAWindowProblem {
public:
    virtual void solve(double threshold, int windowLen,std::vector<approximate_matching::utils::Interval>&) = 0;
};

template<class Integer>
class ThresholdWindowProblem: public  IThresholdAWindowProblem<Integer> {
public:

    ThresholdWindowProblem(AbstractCompleteAMatch<Integer>* aMatch,
                           AbstractPermutation* perm, Integer const *p, int m, Integer const *text, int n, approximate_matching::utils::AbstractScoringScheme* scoringScheme):
                           _aMatch(aMatch),_p(p),_m(m),_text(text),_n(n),_scoringScheme(scoringScheme),permutation(perm){}

    void solve(double threshold, int windowLen,std::vector<approximate_matching::utils::Interval>& result) override {
        if (windowLen > _n) return;
        std::vector<approximate_matching::utils::Interval> innerResult;
        _aMatch->solve(_p,_m,_text,_n,innerResult);
        int m = _m;
        int n = _n;

        auto curCol = windowLen;
        auto curRow = m;

        auto canonicalDecomposition = [this,&m](int i, int j) {
            int rangeSum = 0;
            for (int row = 0; row < permutation->row_size; row++) {
                auto col = permutation->get_col_by_row(row);
                if (row >= i && col < j) rangeSum++;
            }
            return  (j - (i - m) - rangeSum) ;
        };

        auto nextInRow = [this,&m](int i, int j, double value) {
            return  j + 1 - (i - m) - dominance_sum_counting::top_right_arrow::right_move(i, j,j - i + m - value, permutation);
        };
        auto nextInCol = [this,&m](int i, int j, double value) {
            return  j - (i + 1 - m) - dominance_sum_counting::top_right_arrow::down_move(i, j,j - i + m - value, permutation);
        };

        auto curValue = canonicalDecomposition(curRow,curCol);  // O(n)

        auto value  = _scoringScheme->getOriginalScoreFunc(curValue, m, curRow - m , curCol);
//        todo do we need it?
        if (value >= threshold || (innerResult[curCol].start + m > curRow && innerResult[curCol].score >= threshold)) {
            result.emplace_back(curRow - m, curCol, value);
        }

        // goes diagonal via incremental queries
        while (curRow + 1 <= n + m && curCol + 1 <= n) {
            curValue = nextInRow(curRow, curCol, curValue);
            curCol++;
            curValue = nextInCol(curRow, curCol, curValue);
            curRow++;

            auto value  = _scoringScheme->getOriginalScoreFunc(curValue, m, curRow - m , curCol);

            if (value >= threshold || (innerResult[curCol].start + m > curRow && innerResult[curCol].score >= threshold)) {
                result.emplace_back(curRow - m, curCol, value);
            }
        }
    }

private:
    AbstractCompleteAMatch<Integer> *_aMatch;
    AbstractPermutation* permutation;
    Integer const *_p;
    int _m;
    Integer const *_text;
    int _n;
    approximate_matching::utils::AbstractScoringScheme * _scoringScheme;
};

#endif //CPU_WINDOWAMATCH_H
