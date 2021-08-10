

#ifndef CPU_COMPLETEAMATCH_H
#define CPU_COMPLETEAMATCH_H


template<class Integer>
class AbstractCompleteAMatch<Integer> {
        /**
         * For every prefix of text t<0 : j>, the problem asks for the maximum
         * alignment score of pattern p against all possible choices of a suffix t<i : j> from this prefix
         * @return IntArray of positions
         */
        void solve(T*p, int m, T*text,int n,std::vector<Interval> &) = 0;
}

template <class T>
class CompleteAMatchViaSemiLocalTotallyMonotone<T>: public {
    override fun solve(): Array<Pair<Int, Double>> {}
}








#endif //CPU_COMPLETEAMATCH_H
