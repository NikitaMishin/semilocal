

#ifndef CPU_COMPLETEAMATCH_H
#define CPU_COMPLETEAMATCH_H

#include "utils.h"

template<class Integer>
class AbstractCompleteAMatch {
public:
    /**
     * For every prefix of text t<0 : j>, the problem asks for the maximum
     * alignment score of pattern p against all possible choices of a suffix t<i : j> from this prefix
     * @return IntArray of positions
     */
    virtual void solve(Integer const *p, int m, Integer const *text, int n,
                       std::vector<approximate_matching::utils::Interval> &) = 0;
};


template<class Integer>
class SellersCompleteAMatch : public  AbstractCompleteAMatch<Integer> {
public:
    SellersCompleteAMatch(approximate_matching::utils::AbstractScoringScheme *scoringScheme) {
        _scoringScheme = scoringScheme;
    }

    void solve(const Integer *p, int m, const Integer *text, int n,
               std::vector<approximate_matching::utils::Interval> &result) override {
        auto match = _scoringScheme->getMatchScore().toDouble();
        auto mismatch = _scoringScheme->getMismatchScore().toDouble();
        auto gap = _scoringScheme->getGapScore().toDouble();

        auto scoreMatrix = new double[(m + 1) * (n + 1)];
        for (int i = 0; i <  m + 1; ++i) {
            for (int j = 0; j < n + 1; ++j) {
                double fillValue = 0.0;
                if (j == 0) fillValue = i * gap;
                scoreMatrix[i * (n + 1) + j] = fillValue;
            }
        }

        for (int i = 1; i <= m; i++) {
            for (int j = 1; j <= n; ++j) {
                auto cellDiag = scoreMatrix[(i - 1) * (n + 1) + j - 1] + ((p[i - 1] == text[j - 1]) ? match : mismatch);
                auto topCell =  scoreMatrix[(i - 1) * (n + 1) + j] + gap;
                auto leftCell = scoreMatrix[ i * (n + 1) + j - 1] + gap;
                scoreMatrix[i * (n + 1) + j] = std::max(std::max(cellDiag, topCell), leftCell);
            }
        }

        for (int j = 0; j < n; ++j) {
            auto jCap = j;
            auto i = m;
            while (jCap != 0 && i != 0) {
                auto left =     scoreMatrix[ i      * (n+1) + jCap - 1] + gap;
                auto top =      scoreMatrix[(i - 1) * (n+1) + jCap]     + gap;
                auto diagonal = scoreMatrix[(i - 1) * (n+1) + jCap - 1] + ((p[i - 1] == text[jCap - 1]) ? match : mismatch);
                if (top >= left && top >= diagonal) {
                    i--;
                } else if (diagonal >= left && diagonal >= top) {
                    jCap--;
                    i--;
                } else {
                    jCap--;
                }
            }
            result.emplace_back(jCap, j, scoreMatrix[m * (n+1) + j]);
        }
    }

private:
    approximate_matching::utils::AbstractScoringScheme *_scoringScheme;
};

template<class Integer>
class SemiLocalAMatch : public  AbstractCompleteAMatch<Integer> {
public:
    SemiLocalAMatch(approximate_matching::utils::AbstractScoringScheme *scoringScheme) {
        _scoringScheme = scoringScheme;
    }

    void solve(const Integer *p, int m, const Integer *text, int n,
               std::vector<approximate_matching::utils::Interval> &result) override {
        //TODO linear time SMAWK
    }

private:
    approximate_matching::utils::AbstractScoringScheme *_scoringScheme;

    inline double scoreTransformer(double value, int i, int j,int m) {
        _scoringScheme->getOriginalScoreFunc(value,m,i,j);
    }
};

#endif //CPU_COMPLETEAMATCH_H
