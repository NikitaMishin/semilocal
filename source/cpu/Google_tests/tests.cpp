#include <unordered_set>
#include "gtest/gtest.h"
#include "../sequence_generators.h"
#include "../naive_prefix_lcs.h"
#include "../prefix_lcs/bitwise/transposition_network_unbounded_alphabet.h"
#include "../semi_local.h"
#include "../prefix_lcs/bitwise/transposition_network_binary_alphabet.h"
#include "../prefix_lcs/bitwise/encoders_and_decoders.h"

const int MAX_SIZE_B = 1800;
const int ITERATIONS_PER_INPUT = 1;

const int MIN_SIZE_B = 10;
const int STEP = 23;

const int MAX_SIZE_B_BIT = 100;
const int MIN_SIZE_BIT = 1;
const int STEP_BIT = 1;


TEST(prefix_lcs_via_braid_sequential, int_holder) {

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto actual = prefix_lcs_via_braid_sequential<int, int>(seq_a, seq_b);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_sequential, bool_holder) {

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto actual = prefix_lcs_via_braid_sequential<int, bool>(seq_a, seq_b);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_sequential, long_long_holder) {

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto actual = prefix_lcs_via_braid_sequential<int, long long>(seq_a, seq_b);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_mpi_bitwise_operator, long_long_holder_4_thread) {
    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto actual = prefix_lcs_via_braid_mpi_bitwise_operator<int, long long>(seq_a, seq_b, 4);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_mpi_bitwise_operator, char_holder_6_thread) {
    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto actual = prefix_lcs_via_braid_mpi_bitwise_operator<int, char>(seq_a, seq_b, 1);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_mpi_less_operator, long_long_holder_4_thread) {
    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto actual = prefix_lcs_via_braid_mpi_less_operator<int, long long>(seq_a, seq_b, 4);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_mpi_less_operator, char_holder_6_thread) {
    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto actual = prefix_lcs_via_braid_mpi_less_operator<int, char>(seq_a, seq_b, 1);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(sticky_braid_mpi, int_holder_4threads) {

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = sticky_braid_sequential<int, int>(seq_a, seq_b);
                auto actual = sticky_braid_mpi<int, int>(seq_a, seq_b, 4);
                for (int j = 0; j < seq_a.size() + seq_b.size(); ++j) {
                    EXPECT_EQ(expected[i], actual[i]);
                }

            }
        }
    }
}

TEST(sticky_braid_mpi, long_long_holder_6threads) {

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 26);
                auto seq_b = gen_vector_seq(b_size, 26);
                auto expected = sticky_braid_sequential<int, long long>(seq_a, seq_b);
                auto actual = sticky_braid_mpi<int, long long>(seq_a, seq_b, 6);
                for (int j = 0; j < seq_a.size() + seq_b.size(); ++j) {
                    EXPECT_EQ(expected[i], actual[i]);
                }

            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_binary, unsigned_int_holder) {
    typedef unsigned int wordType;
    auto size = sizeof(unsigned int) * 8;

    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 2);
                auto seq_b = gen_vector_seq(b_size, 2);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_binary(a.first.first, a.first.second, a.second,
                                                               b.first.first, b.first.second, b.second);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_binary, size_t_holder) {
    typedef size_t wordType;
    auto size = sizeof(size_t) * 8;

    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 2);
                auto seq_b = gen_vector_seq(b_size, 2);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_binary(a.first.first, a.first.second, a.second,
                                                               b.first.first, b.first.second, b.second);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_binary, unsigned_char_holder) {
    typedef unsigned char wordType;
    auto size = sizeof(unsigned char) * 8;

    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 2);
                auto seq_b = gen_vector_seq(b_size, 2);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_binary(a.first.first, a.first.second, a.second,
                                                               b.first.first, b.first.second, b.second);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_binary_mpi, unsigned_int_holder_4threads) {
    typedef unsigned int wordType;
    auto size = sizeof(unsigned int) * 8;

    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 2);
                auto seq_b = gen_vector_seq(b_size, 2);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_binary_mpi(a.first.first, a.first.second, a.second,
                                                                   b.first.first, b.first.second, b.second, 4);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_binary_mpi, size_t_holder_4threads) {
    typedef size_t wordType;
    auto size = sizeof(size_t) * 8;

    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size; b_size += STEP_BIT * size) {
        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 2);
                auto seq_b = gen_vector_seq(b_size, 2);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_binary_mpi(a.first.first, a.first.second, a.second,
                                                                   b.first.first, b.first.second, b.second, 4);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_binary_mpi, unsigned_char_holder_4threads) {
    typedef unsigned char wordType;
    auto size = sizeof(unsigned char) * 8;

    for (int b_size = MIN_SIZE_BIT * size; b_size < MAX_SIZE_B_BIT * size * 2; b_size += STEP_BIT * size) {
        for (int a_size = size; a_size <= b_size; a_size += STEP_BIT * size) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 2);
                auto seq_b = gen_vector_seq(b_size, 2);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_binary_mpi(a.first.first, a.first.second, a.second,
                                                                   b.first.first, b.first.second, b.second, 4);
                EXPECT_EQ(expected, actual);

            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_4symbol, unsigned_int_holder) {
    typedef unsigned int wordType;

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 4);
                auto seq_b = gen_vector_seq(b_size, 4);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_4symbol(a.first.first, a.first.second, a.second,
                                                                b.first.first, b.first.second, b.second);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_4symbol, size_tholder) {
    typedef size_t wordType;

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 4);
                auto seq_b = gen_vector_seq(b_size, 4);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_4symbol(a.first.first, a.first.second, a.second,
                                                                b.first.first, b.first.second, b.second);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_bits_4symbol, unsigned_char_holder) {
    typedef unsigned char wordType;

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 4);
                auto seq_b = gen_vector_seq(b_size, 4);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_4symbol<wordType>(a.first.first, a.first.second, a.second,
                                                                          b.first.first, b.first.second, b.second);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_bits_4symbol_mpi, unsigned_int_holder_4threads) {
    typedef unsigned int wordType;

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 4);
                auto seq_b = gen_vector_seq(b_size, 4);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_4symbol_mpi(a.first.first, a.first.second, a.second,
                                                                    b.first.first, b.first.second, b.second, 4);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}


TEST(prefix_lcs_via_braid_bits_4symbol_mpi, size_tholder_4threads) {
    typedef size_t wordType;

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 4);
                auto seq_b = gen_vector_seq(b_size, 4);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_4symbol_mpi(a.first.first, a.first.second, a.second,
                                                                    b.first.first, b.first.second, b.second, 4);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}

TEST(prefix_lcs_via_braid_bits_4symbol_mpi, unsigned_char_holder_4threads) {
    typedef unsigned char wordType;

    for (int b_size = MIN_SIZE_B; b_size < MAX_SIZE_B; b_size += STEP) {
        for (int a_size = 1; a_size <= b_size; a_size += STEP) {
            for (int i = 0; i < ITERATIONS_PER_INPUT; ++i) {
                auto seq_a = gen_vector_seq(a_size, 4);
                auto seq_b = gen_vector_seq(b_size, 4);
                auto expected = prefix_lcs_sequential(seq_a, seq_b);
                auto mappers = encode_alphabet<int, wordType>(std::unordered_set<int>({2, 3, 1, 0}));
                auto a = encode_reverse<int, wordType>(seq_a, &mappers.first, &mappers.second);
                auto b = encode<int, wordType>(seq_b, &mappers.first, &mappers.second);
                auto actual = prefix_lcs_via_braid_bits_4symbol_mpi<wordType>(a.first.first, a.first.second, a.second,
                                                                              b.first.first, b.first.second, b.second,
                                                                              4);
                EXPECT_EQ(expected, actual);
            }
        }
    }
}