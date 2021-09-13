
#ifndef CPU_THRESHOLDAMATCH_H
#define CPU_THRESHOLDAMATCH_H

#include "utils.h"

namespace approximate_matching {


    template<class Integer>
    class AbstractThresholdAMatch {
    public:
        /**
         * For a given threshold in [0,1] according to scoring scheme find
         */
        virtual void
        find(double threshold, Integer const *p, int m, Integer const *text, int n, std::vector<utils::Interval> &) = 0;
    };

    template<class Integer>
    class ThresholdAMatch {
    public:
        ThresholdAMatch(AbstractCompleteAMatch<Integer> *aMatch) {
            _aMatch = aMatch;
        }

        void find(double threshold, Integer const *p, int m, Integer const *text, int n,
                  std::vector<utils::Interval> &result) {
            std::vector<utils::Interval> innerResult;
            _aMatch->solve(p, m, text, n, innerResult);

            auto j = innerResult.size() - 1;
            //will go backwards with filtering
            std::reverse(result.begin(), result.end());
            while (j > 0) {
                if (innerResult[j].score >= threshold) {
                    result.push_back(innerResult[j]);
                    break;
                }
                j--;
            }

            //actual filtration
            while (j > 0) {
                // or remove if clase todo what is better heuristic
                auto &last = result[result.size() - 1];
                if (innerResult[j].score >= threshold &&
                    (innerResult[j].end > last.start && innerResult[j].score >= last.score) &&
                    innerResult[j].start >= last.start) {
                    last.start = innerResult[j].start;
                    last.end = innerResult[j].end;
                    last.score = innerResult[j].score;
                } else if (innerResult[j].score >= threshold && j <= result[result.size() - 1].start)
                    result.push_back(innerResult[j]);
                j--;
            }

            std::reverse(result.begin(), result.end());
        }

    private:
        AbstractCompleteAMatch<Integer> *_aMatch;
    };
}


#endif //CPU_THRESHOLDAMATCH_H
