import os
from argparse import ArgumentParser
from os import listdir

from utils.logger import Logger
from utils.runners import CombingRunner
from utils.testing import run_tests, compile_programs
from utils.tests import CombingTest

# CXX_COMPILER_PATH = '/usr/bin/g++-10'
CXX_COMPILER_PATH = '/usr/bin/g++'

SINGLE_THREADED_SOLUTIONS = [

]

MULTI_THREADED_SOLUTIONS = [
    'prefix_lcs_semi_bit_parallel_1form',
    'prefix_lcs_semi_bit_parallel_2form',
    'prefix_lcs_semi_bit_parallel_old',
]

SOLUTIONS_FOLDER = 'combing_solutions'  # where we put our ready to run implementations

default_logger = Logger('fig8910_log.txt')
build_logger = Logger('fig8910_build_log.txt')


def boolean_string(s):
    if s not in {'False', 'True'}:
        raise ValueError('Not a valid boolean string')
    return s == 'True'


if __name__ == '__main__':
    arg_parser = ArgumentParser()
    arg_parser.add_argument('max_thds', type=int, help='Threads')
    arg_parser.add_argument('tests', type=str, help='Path to dataset with *.fna files, our specific format')

    arg_parser.add_argument('resultcsvfile', type=str, help='csv file')
    arg_parser.add_argument('repeats', type=int, help='repeats')

    args = arg_parser.parse_args()

    runners = compile_programs(SINGLE_THREADED_SOLUTIONS, MULTI_THREADED_SOLUTIONS, 1, args.max_thds + 1, CombingRunner,
                               build_logger, SOLUTIONS_FOLDER, CXX_COMPILER_PATH)

    test_cases = []
    for name in listdir(args.tests):
        with open(os.path.join(args.tests, name), 'r') as f:
            _, size = f.readline(), int(f.readline())
            test_cases.append((size, os.path.join(args.tests, name)))


    tests_pairs = [(100000,100000), (200000,200000), (400000,400000),(600000,600000),(800000,800000),(1000000,1000000)]

    tests = []

    for i in range(len(test_cases)):
        for j in range(len(test_cases)):
            x, y = test_cases[i], test_cases[j]
            if x[1] != y[1] and (x[0], y[0]) in tests_pairs:
                tests.append((x[0], CombingTest(x[1], y[1])))

    tests.sort(key=lambda x: x[0])

    _, tests = zip(*tests)
    tests = list(tests)

    run_tests(runners, tests, args.repeats, args.resultcsvfile, default_logger, True)
