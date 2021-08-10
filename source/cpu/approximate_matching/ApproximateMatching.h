
#ifndef CPU_APPROXIMATEMATCHING_H
#define CPU_APPROXIMATEMATCHING_H

namespace approximate_matching {

    struct Interval {
        double score;
        int start;
        int end;
    };

    template<class Integer>
    class AbstractApproximateMatching {
    public:
        /**
         * For a given threshold in [0,1] according to scoring scheme find
         */
        virtual void find(T *p, int m, T *text, int n, std::vector <Interval> &) = 0;
    }
}



#endif //CPU_APPROXIMATEMATCHING_H
