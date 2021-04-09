import os
from argparse import ArgumentParser
from os import listdir

from utils.logger import Logger
from utils.runners import BraidMulRunner, CombingRunner
from utils.testing import run_tests, compile_programs
from utils.tests import CombingTest

# CXX_COMPILER_PATH = '/usr/bin/g++-10'
CXX_COMPILER_PATH = '/usr/bin/g++'

SINGLE_THREADED_SOLUTIONS = [
    'semi_local_naive_iterative',
    'prefix_lcs',
    'prefix_lcs_skewed',

    'semi_local_parallel_withoutif_iterative',
    'semi_local_parallel_withif_iterative',
    'semi_local_parallel_1and3_combined_iterative'
]

MULTI_THREADED_SOLUTIONS = [
    # 'semi_local_parallel_iterative',
]

SOLUTIONS_FOLDER = 'combing_solutions'  # where we put our ready to run implementations


default_logger = Logger('fig4_log.txt')
build_logger = Logger('fig4_build_log.txt')


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'



if __name__ == '__main__':
    arg_parser = ArgumentParser()
    arg_parser.add_argument('tests', type=str, help='Path to dataset with *.fna files, our specific format')
    arg_parser.add_argument('resultcsvfile', type=str, help='csv file')
    arg_parser.add_argument('repeats', type=int, help='repeats')

    arg_parser.add_argument('real_data', type=boolean_string,
                            help='tests is real data')

    args = arg_parser.parse_args()

    runners = compile_programs(SINGLE_THREADED_SOLUTIONS, MULTI_THREADED_SOLUTIONS, 1, 1, CombingRunner,
                               build_logger, SOLUTIONS_FOLDER, CXX_COMPILER_PATH)

    test_cases = []
    for name in listdir(args.tests):
        with open(os.path.join(args.tests, name), 'r') as f:
            _, size = f.readline(), int(f.readline())
            test_cases.append((size, os.path.join(args.tests, name)))

    pair_tests = [(7742, 9663), (13517, 16945), (33356, 39245), (43548, 50884), (57623, 59815), (74437, 91721),
                  (124884, 134226)]
    tests = []

    if args.real_data:
        for i in range(len(test_cases)):
            for j in range(len(test_cases)):
                x, y = test_cases[i], test_cases[j]
                if x[1] != y[1] and (x[0], y[0]) in pair_tests:
                    print((x[0], y[0]))
                    tests.append((x[0], CombingTest(x[1], y[1])))
    else:
        for pair_1 in test_cases:
            test_case_1 = pair_1[1]
            size_test_case_1 = pair_1[0]
            tests.extend([(size_test_case_1, CombingTest(test_case_1, test_name)) for size, test_name in test_cases if
                          size == size_test_case_1 and test_case_1 != test_name])

    tests.sort(key=lambda x: x[0])


    _, tests = zip(*tests)
    tests = list(tests)

    run_tests(runners, tests, args.repeats, args.resultcsvfile, default_logger, True)
