#pragma once


template<class Input>
class LCSStrategy {
public:
    virtual int compute(const Input *a, int sizeA, const Input *b, int sizeB);
};

template<class Input>
class DummyLCS : public LCSStrategy<Input> {
public:
    int compute(const Input *a, int sizeA, const Input *b, int sizeB) override {
        int arr[sizeA + 1][sizeB + 1];
        auto m = sizeA + 1;
        auto n = sizeB + 1;

        for (auto j = 0; j < n; j++) arr[0][j] = 0;
        for (auto i = 0; i < m; i++) arr[i][0] = 0;

        for (int i = 1; i < m; ++i) {
            for (int j = 1; j < n; ++j) {
                arr[i][j] = std::max(std::max(arr[i - 1][j], arr[i][j - 1]), (a[i - 1] == b[j - 1]) ? arr[i - 1][j - 1] + 1 : arr[i - 1][j - 1]);
            }
        }
        return arr[m - 1][n - 1];
    }


};