//
// Created by garrancha on 10.08.2021.
//

#ifndef CPU_UTILS_H
#define CPU_UTILS_H

#include <string>
#include <string_view>
#include <cmath>
#include "../unit_monge_mult/matrices.h"
#include "../unit_monge_mult/dominance_sum_queries.h"

#define SPECIAL_SYMBOL -1
namespace approximate_matching {
    namespace utils {


        template<class Integer>
        void blownup(Integer *string, int string_size, int v, int mu, Integer *result) {
            auto index = 0;
            for (int k = 0; k < string_size; k++) {
                auto character = string[k];
                for (int i = 0; i < mu; i++, index++) result[index] = SPECIAL_SYMBOL;
                for (int i = 0; i < v - mu; i++, index++) result[index] = character;
            }
        }

        int gcd(int g, int e) {
            auto a = abs(g);
            auto b = abs(e);
            if (a == 0 || b == 0) return std::max(a, b);
            while (b != 0) {
                int t = b;
                b = a % b;
                a = t;
            }
            return a;
        }

        class Fraction {
        private:
            int _numerator;
            int _denominator;
            int sign;
        public:
            Fraction(int numerator, int denominator) {
                if (denominator == 0) throw std::runtime_error("Denominator mustn't equals to zero");

                int gcd_value = gcd(numerator, denominator);
                _numerator = abs(numerator / gcd_value);
                _denominator = abs(denominator / gcd_value);
                auto sign_a = numerator > 0;
                auto sign_b = denominator > 0;
                if (sign_a == sign_b) {
                    sign = 1;
                } else if (numerator == 0 || denominator == 0) {
                    sign = 0;
                } else {
                    sign = -1;
                }
            };

            int getNumerator() {
                return _numerator;
            }

            int getDenominator() {
                return _denominator;
            }

            double toDouble() const {
                return sign * double(_numerator) / _denominator;
            }

            bool operator==(Fraction &other) const {
                return this->sign == other.sign && other._numerator == _numerator && other._denominator == _denominator;
            }

            template<class Word>
            Word hashCode() {
                Word hash = 7;
                hash = 31 * hash + _numerator;
                hash = 31 * hash + _denominator;
                hash = 31 * hash + sign;
                return hash;
            }

            std::string toString() const {
                switch (sign) {
                    case 1:
                        return std::to_string(_numerator) + std::to_string(_denominator);
                    case -1:
                        return "-" + std::to_string(_numerator) + std::to_string(_denominator);
                    default:
                        return "0";
                }
            }

            Fraction operator+(Fraction &other) const {
                return {
                        sign * other._denominator * _numerator + other.sign * other._numerator * _denominator,
                        _denominator * other._denominator
                };
            }

            Fraction operator*(Fraction &other) const {
                if (_numerator == 0 || other._numerator == 0) return {0, 1};
                return {
                        sign * other.sign * _numerator * other._numerator,
                        _denominator * other._denominator
                };
            }

            Fraction operator-(Fraction &other) const {
                return {
                        sign * other._denominator * _numerator - other.sign * other._numerator * _denominator,
                        _denominator * other._denominator
                };
            }

            Fraction operator-(Fraction other) const {
                return {
                        sign * other._denominator * _numerator - other.sign * other._numerator * _denominator,
                        _denominator * other._denominator
                };
            }

            Fraction operator/(Fraction &other) const {
                if (other._numerator == 0) throw std::runtime_error("Division by zero");
                if (_numerator == 0 || other._denominator == 0) return {0, 1};
                return {sign * other.sign * _numerator * other._denominator,
                        _denominator * other._numerator
                };
            }

            Fraction operator/(Fraction other) const {
                if (other._numerator == 0) throw std::runtime_error("Division by zero");
                if (_numerator == 0 || other._denominator == 0) return {0, 1};
                return {sign * other.sign * _numerator * other._denominator,
                        _denominator * other._numerator
                };
            }

            Fraction operator+(int other) {
                auto frac = Fraction(other, 1);
                return (*this) + frac;
            }

            Fraction operator*(int other) {
                auto frac = Fraction(other, 1);
                return (*this) * frac;
            }


        };

        class AbstractScoringScheme {
        public:
            AbstractScoringScheme() = default;;

            virtual Fraction getNormalizedMatchScore() = 0;

            virtual Fraction getNormalizedMismatchScore() = 0;

            virtual Fraction getNormalizedGapScore() = 0;

            virtual Fraction getMatchScore() = 0;

            virtual Fraction getMismatchScore() = 0;

            virtual Fraction getGapScore() = 0;

            virtual double getOriginalScoreFunc(double value, int m, int i, int j) = 0;

            ~AbstractScoringScheme() = default;
        };

        class RegularScoringScheme : public AbstractScoringScheme {
        public:
            explicit RegularScoringScheme(Fraction frac) : mismatch(frac), AbstractScoringScheme() {
                if (frac.toDouble() > 1.0 || frac.toDouble() < -1.0) std::runtime_error("Could not be higher then 1");
            };

            Fraction getNormalizedMatchScore() override { return match; };

            Fraction getNormalizedMismatchScore() override { return mismatch; }

            Fraction getNormalizedGapScore() override { return gap; }

            Fraction getMatchScore() override { return match; }

            Fraction getMismatchScore() override { return mismatch; }

            Fraction getGapScore() override { return gap; }

            double getOriginalScoreFunc(double value, int m, int i, int j) override { return value; }

        private:
            Fraction match = Fraction(1, 1);
            Fraction mismatch;
            Fraction gap = Fraction(0, 1);


        };


        class FixedScoringScheme : public AbstractScoringScheme {
        public:
            FixedScoringScheme(Fraction match, Fraction mismatch, Fraction gap) : _match(match), _mismatch(mismatch),
                                                                                  _gap(gap),
                                                                                  normalizedMismatch(
                                                                                          (mismatch - (gap * 2)) /
                                                                                          (match - (gap * 2))) {}

            Fraction getNormalizedMatchScore() override { return {0, 1}; }

            Fraction getNormalizedMismatchScore() override { return normalizedMismatch; }

            Fraction getNormalizedGapScore() override { return {0, 1}; }


            Fraction getMatchScore() override { return _match; };

            Fraction getMismatchScore() override { return _mismatch; }

            Fraction getGapScore() override { return _gap; }

            double getOriginalScoreFunc(double value, int m, int i, int j) override {
                return value * (_match - (_gap * 2)).toDouble() + ((m + j - i)) * _gap.toDouble();
            }

            double getReverseScore(double originalValue, int m, int i, int j) {
                return (originalValue - (m + j - i) * _gap.toDouble()) / ((_match - (_gap * 2)).toDouble());
            }

        private:
            Fraction _match;
            Fraction _mismatch;
            Fraction _gap;
            Fraction normalizedMismatch;
        };


        struct Interval {
            int start;
            int end;
            double score;

            Interval() {};

            Interval(int start_, int end_, int score_) : start(start_), end(end_), score(score_) {}

            inline int len() const {
                return end - start;
            }

            bool operator==(const Interval &other) const {
                return other.start == start && other.end == end && other.score == score;
            }

            friend auto operator<<(std::ostream &os, Interval const &m) -> std::ostream & {
                return os << m.score << "," << m.start << ":" << m.end;
            }
        };


        class SemiLocalStringSubstringWrapper {
        public:

            int p_size;
            int v;
            Permutation *perm;
            AbstractScoringScheme *scheme;

            SemiLocalStringSubstringWrapper(Permutation *permutation, AbstractScoringScheme *scheme_,
                                            int pattern_size, int v_value) : perm(permutation), p_size(pattern_size),
                                                                             v(v_value), scheme(scheme_) {}


            int dominanceSumRandom(int i, int j) {
                i += p_size;
                i *= v;
                j *= v;
                int rangeSum = 0;
                for (int row = 0; row < perm->row_size; row++) {
                    auto col = perm->get_col_by_row(row);
                    if (row >= i && col < j) {
                        rangeSum++;
                    }
                }
                return rangeSum;
            };

            int dominanceSumForMoveUp(int i, int j, int value) {
                i += p_size;
                i *= v;
                j *= v;

                for (int l = 0; l < v; ++l) {
                    value = dominance_sum_counting::top_right_arrow::up_move(i, j, value, perm);
                    i--;
                }
                return value;
            };

            int dominanceSumForMoveDown(int i, int j, int value) {
                i += p_size;
                i *= v;
                j *= v;

                for (int l = 0; l < v; ++l) {
                    value = dominance_sum_counting::top_right_arrow::down_move(i, j, value, perm);
                    i++;
                }

                return value;
            };

            int dominanceSumForMoveLeft(int i, int j, int value) {
                i += p_size;
                i *= v;
                j *= v;
                for (int l = 0; l < v; ++l) {
                    value = dominance_sum_counting::top_right_arrow::left_move(i, j, value, perm);
                    j--;
                }
                return value;
            };

            double canonicalDecompositionWithKnown(int i, int j, int rangeSum) {
                i += p_size;
                return j - (i - p_size) - double(rangeSum) / v;
            };
            // 5 - 2.5 =2.5
            //


            double originalScore(int i, int j, int dominanceSum) {
                return scheme->getOriginalScoreFunc(canonicalDecompositionWithKnown(i, j, dominanceSum), p_size,
                                                    i, j);
            };

            int getRangeSumFromReverse(int i, int j, double result) {
                // original flow:
                // t =  j - i - rangeSum / v
                // result = t* ( a - 2b) + (m+j-i)b
                // backward flow:
                // rangeSum = v(j - i - t)
                // t = (result - (m + j - i)b ) / (a - 2b)
                // substitution
                // rangeSum = v (j-i - (result - (m+j-i)b)/(a-2b)  )
                auto a = scheme->getMatchScore().toDouble();
                auto b = scheme->getGapScore().toDouble();
                return std::round(v * (j - i - (result - (p_size + j - i) * b) / (a - 2 * b)));
            }


        };

    }


}


#endif //CPU_UTILS_H
