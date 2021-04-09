#include <unordered_set>
#include "gtest/gtest.h"
#include "../sequence_generators.h"
#include "../naive_prefix_lcs.h"
#include "../prefix_lcs/bitwise/transposition_network_unbounded_alphabet.h"
#include "../semi_local.h"
#include "../prefix_lcs/bitwise/transposition_network_binary_alphabet.h"
#include "../prefix_lcs/bitwise/encoders_and_decoders.h"

const int MAX_SIZE_B = 100 * 3;
const int ITERATIONS_PER_INPUT = 1;

const int MIN_SIZE_B = 10;
const int STEP = 1;

const int MAX_SIZE_B_BIT = 100;
const int MIN_SIZE_BIT = 1;
const int STEP_BIT = 1;


/**
 * Testing correctness of skewed version of prefix lcs
 */
TEST(prefix_lcs, prefix_lcs_vs_skewed_prefix) {

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                for (int alphabet = 1; alphabet < 25; ++alphabet) {
                    auto seq_a = gen_vector_seq<int>(a_size, alphabet);
                    auto seq_b = gen_vector_seq<int>(b_size, alphabet);
                    auto a = new int[a_size];
                    auto b = new int[b_size];
                    for (int j = 0; j < a_size; ++j) a[j] = seq_a[j];
                    for (int j = 0; j < b_size; ++j) b[j] = seq_b[j];


                    auto expected = prefix_lcs_sequential(a, a_size, b, b_size);
                    auto actual = prefix_lcs_sequential_skewed(a, a_size, b, b_size);

                    delete[] a;
                    delete[] b;

                    EXPECT_EQ(expected, actual);
                }

            }
        }
    }
}


TEST(semi_local_lcs, antidiagonal_vs_naive) {

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq<int>(a_size, 26);
                auto seq_b = gen_vector_seq<int>(b_size, 26);

                auto a = new int[a_size];
                auto b = new int[b_size];
                for (int j = 0; j < a_size; ++j) a[j] = seq_a[j];
                for (int j = 0; j < b_size; ++j) b[j] = seq_b[j];

                auto should = Permutation(a_size + b_size, a_size + b_size);
                auto actual = Permutation(a_size + b_size, a_size + b_size);

                semi_local::sticky_braid_sequential<int, false>(should, a, a_size, b, b_size);
                semi_local::sticky_braid_mpi<int, false, true>(actual, a, a_size, b, b_size, 1);


                EXPECT_EQ(1, should.is_equal_to(actual));
            }
        }
    }
}

TEST(semi_local_lcs, down_top_approach) {
    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
    precalc(map, 5);


    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP + 11) {
        for (int a_size = 10; a_size <= b_size; a_size += STEP + 11) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                for (int t = 1; t < 10; ++t) {

                    for (int k = 1; k < 10; ++k) {
                        auto seq_a = gen_vector_seq<int>(a_size, 26);
                        auto seq_b = gen_vector_seq<int>(b_size, 26);

                        auto a = new int[a_size];
                        auto b = new int[b_size];
                        for (int j = 0; j < a_size; ++j) a[j] = seq_a[j];
                        for (int j = 0; j < b_size; ++j) b[j] = seq_b[j];

                        auto should = Permutation(a_size + b_size, a_size + b_size);
                        auto actual = Permutation(a_size + b_size, a_size + b_size);

                        semi_local::sticky_braid_sequential<int, false>(should, a, a_size, b, b_size);
                        semi_local::semi_local_down_to_top<int, false>(actual, a, a_size, b, b_size, map, t, k);


                        EXPECT_EQ(1, should.is_equal_to(actual));

                    }
                }

            }
        }
    }
}


TEST(semi_local_lcs, hybrid) {
    auto map = std::unordered_map<int, std::unordered_map<long long, std::unordered_map<long long, std::vector<std::pair<int, int>>>>>();
    precalc(map, 5);


    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 10; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq<int>(a_size, 26);
                auto seq_b = gen_vector_seq<int>(b_size, 26);

                auto a = new int[a_size];
                auto b = new int[b_size];
                for (int j = 0; j < a_size; ++j) a[j] = seq_a[j];
                for (int j = 0; j < b_size; ++j) b[j] = seq_b[j];

                auto should = Permutation(a_size + b_size, a_size + b_size);
                auto actual = Permutation(a_size + b_size, a_size + b_size);

                semi_local::sticky_braid_sequential<int, false>(should, a, a_size, b, b_size);
                semi_local::hybrid<int, false, true, true>(actual, a, a_size, b, b_size, map, 1, 0, 0, 7, 0);


                EXPECT_EQ(1, should.is_equal_to(actual));

            }
        }

    }
}

TEST(llcs_2symbol_smart_combing, llcs_2symbol_smart_combing_smart) {

    for (int b_size = 32; b_size < 2000; b_size += 32) {
        for (int a_size = 32; a_size <= b_size; a_size += 32) {
            for (int i = 0; i < 2; ++i) {
                auto seq_a = gen_vector_seq<int>(a_size, 2);
                auto seq_b = gen_vector_seq<int>(b_size, 2);


                auto mappers = encode_alphabet<int,unsigned int>(std::unordered_set<int>({0,1}));
                auto a = encode_reverse<int, unsigned int>(seq_a, &mappers.first, &mappers.second,2);
                auto b = encode<int, unsigned int>(seq_b,&mappers.first,&mappers.second,2);

                auto a_int = new int[a_size];
                auto b_int = new int[b_size];
                for (int j = 0; j < a_size; ++j) a_int[j] = seq_a[j];
                for (int j = 0; j < b_size; ++j) b_int[j] = seq_b[j];

                auto expected = prefix_lcs_sequential(a_int, a_size, b_int, b_size);

                auto actual = prefix_lcs_via_semi_local::binary::llcs_2symbol_smart_combing<unsigned int>(
                        a.first.first,a.first.second, b.first.first,b.first.second, a.second,1, i % 2);

                EXPECT_EQ(expected,actual);

            }
        }

    }
}

TEST(llcs_2symbol_smart_combing, llcs_2symbol_smart_combing_naive) {


    for (int b_size = 32; b_size < 2000; b_size += 32) {
        for (int a_size = 32; a_size <= b_size; a_size += 32) {
            for (int i = 0; i < 2; ++i) {
                auto seq_a = gen_vector_seq<int>(a_size, 2);
                auto seq_b = gen_vector_seq<int>(b_size, 2);


                auto mappers = encode_alphabet<int,unsigned int>(std::unordered_set<int>({0,1}));
                auto a = encode_reverse<int, unsigned int>(seq_a, &mappers.first, &mappers.second,2);
                auto b = encode<int, unsigned int>(seq_b,&mappers.first,&mappers.second,2);

                auto a_int = new int[a_size];
                auto b_int = new int[b_size];
                for (int j = 0; j < a_size; ++j) a_int[j] = seq_a[j];
                for (int j = 0; j < b_size; ++j) b_int[j] = seq_b[j];

                auto expected = prefix_lcs_sequential(a_int, a_size, b_int, b_size);

                auto actual = prefix_lcs_via_semi_local::binary::llcs_2symbol_naive_combing<unsigned int>(
                        a.first.first,a.first.second, b.first.first,b.first.second, a.second,1);

                EXPECT_EQ(expected,actual);

            }
        }

    }
}


//
//TEST(prefix_lcs_via_braid_sequential, int_holder) {
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto actual = prefix_lcs_via_braid_sequential<int, int>(seq_a, seq_b);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//TEST(prefix_lcs_via_braid_sequential, bool_holder) {
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto actual = prefix_lcs_via_braid_sequential<int, bool>(seq_a, seq_b);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//

//
//TEST(prefix_lcs_via_braid_mpi_bitwise_operator, long_long_holder_4_thread) {
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto actual = prefix_lcs_via_braid_mpi_bitwise_operator<int, long long>(seq_a, seq_b, 4);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//TEST(prefix_lcs_via_braid_mpi_bitwise_operator, char_holder_6_thread) {
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto actual = prefix_lcs_via_braid_mpi_bitwise_operator<int, char>(seq_a, seq_b, 1);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_mpi_less_operator, long_long_holder_4_thread) {
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto actual = prefix_lcs_via_braid_mpi_less_operator<int, long long>(seq_a, seq_b, 4);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//TEST(prefix_lcs_via_braid_mpi_less_operator, char_holder_6_thread) {
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto actual = prefix_lcs_via_braid_mpi_less_operator<int, char>(seq_a, seq_b, 1);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(sticky_braid_mpi, int_holder_4threads) {
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = sticky_braid_sequential<int, int>(seq_a, seq_b);
//                auto actual = sticky_braid_mpi<int, int>(seq_a, seq_b, 4);
//                for (int j = 0; j < seq_a.size() + seq_b.size(); ++j) {
//                    EXPECT_EQ(expected[i], actual[i]);
//                }
//
//            }
//        }
//    }
//}
//
//TEST(sticky_braid_mpi, long_long_holder_6threads) {
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 26);
//                auto seq_b = gen_vector_seq(b_size, 26);
//                auto expected = sticky_braid_sequential<int, long long>(seq_a, seq_b);
//                auto actual = sticky_braid_mpi<int, long long>(seq_a, seq_b, 6);
//                for (int j = 0; j < seq_a.size() + seq_b.size(); ++j) {
//                    EXPECT_EQ(expected[i], actual[i]);
//                }
//
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_binary, unsigned_int_holder) {
//    typedef unsigned int wordType;
//    auto size = sizeof(unsigned int) * 8;
//
//    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
//        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 2);
//                auto seq_b = gen_vector_seq(b_size, 2);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_binary(a.first.first, a.first.second, a.second,
//                                                               b.first.first, b.first.second, b.second);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_binary, size_t_holder) {
//    typedef size_t wordType;
//    auto size = sizeof(size_t) * 8;
//
//    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
//        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 2);
//                auto seq_b = gen_vector_seq(b_size, 2);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_binary(a.first.first, a.first.second, a.second,
//                                                               b.first.first, b.first.second, b.second);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_binary, unsigned_char_holder) {
//    typedef unsigned char wordType;
//    auto size = sizeof(unsigned char) * 8;
//
//    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
//        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 2);
//                auto seq_b = gen_vector_seq(b_size, 2);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_binary(a.first.first, a.first.second, a.second,
//                                                               b.first.first, b.first.second, b.second);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_binary_mpi, unsigned_int_holder_4threads) {
//    typedef unsigned int wordType;
//    auto size = sizeof(unsigned int) * 8;
//
//    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
//        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 2);
//                auto seq_b = gen_vector_seq(b_size, 2);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_binary_mpi(a.first.first, a.first.second, a.second,
//                                                                   b.first.first, b.first.second, b.second, 4);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_binary_mpi, size_t_holder_4threads) {
//    typedef size_t wordType;
//    auto size = sizeof(size_t) * 8;
//
//    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
//        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 2);
//                auto seq_b = gen_vector_seq(b_size, 2);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_binary_mpi(a.first.first, a.first.second, a.second,
//                                                                   b.first.first, b.first.second, b.second, 4);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_binary_mpi, unsigned_char_holder_4threads) {
//    typedef unsigned char wordType;
//    auto size = sizeof(unsigned char) * 8;
//
//    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size * 2; b_size += STEP_BIT * size) {
//        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 2);
//                auto seq_b = gen_vector_seq(b_size, 2);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_binary_mpi(a.first.first, a.first.second, a.second,
//                                                                   b.first.first, b.first.second, b.second, 4);
//                EXPECT_EQ(expected, actual);
//
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_4symbol, unsigned_int_holder) {
//    typedef unsigned int wordType;
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 4);
//                auto seq_b = gen_vector_seq(b_size, 4);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_4symbol(a.first.first, a.first.second, a.second,
//                                                                b.first.first, b.first.second, b.second);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_4symbol, size_tholder) {
//    typedef size_t wordType;
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 4);
//                auto seq_b = gen_vector_seq(b_size, 4);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_4symbol(a.first.first, a.first.second, a.second,
//                                                                b.first.first, b.first.second, b.second);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//TEST(prefix_lcs_via_braid_bits_4symbol, unsigned_char_holder) {
//    typedef unsigned char wordType;
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 4);
//                auto seq_b = gen_vector_seq(b_size, 4);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_4symbol<wordType>(a.first.first, a.first.second, a.second,
//                                                                          b.first.first, b.first.second, b.second);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//TEST(prefix_lcs_via_braid_bits_4symbol_mpi, unsigned_int_holder_4threads) {
//    typedef unsigned int wordType;
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 4);
//                auto seq_b = gen_vector_seq(b_size, 4);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_4symbol_mpi(a.first.first, a.first.second, a.second,
//                                                                    b.first.first, b.first.second, b.second, 4);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//
//TEST(prefix_lcs_via_braid_bits_4symbol_mpi, size_tholder_4threads) {
//    typedef size_t wordType;
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 4);
//                auto seq_b = gen_vector_seq(b_size, 4);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_4symbol_mpi(a.first.first, a.first.second, a.second,
//                                                                    b.first.first, b.first.second, b.second, 4);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}
//
//TEST(prefix_lcs_via_braid_bits_4symbol_mpi, unsigned_char_holder_4threads) {
//    typedef unsigned char wordType;
//
//    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
//        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
//            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
//                auto seq_a = gen_vector_seq(a_size, 4);
//                auto seq_b = gen_vector_seq(b_size, 4);
//                auto expected = prefix_lcs_sequential(seq_a, seq_b);
//                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
//                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
//                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
//                auto actual = prefix_lcs_via_braid_bits_4symbol_mpi<wordType>(a.first.first, a.first.second, a.second,
//                                                                              b.first.first, b.first.second, b.second,
//                                                                              4);
//                EXPECT_EQ(expected, actual);
//            }
//        }
//    }
//}