from argparse import ArgumentParser

from utils.logger import Logger
from utils.runners import BraidMulRunner
from utils.testing import run_tests, compile_programs
from utils.tests import BraidMulTest

# /usr/bin/g++-10
CXX_COMPILER_PATH = '/usr/bin/g++-10'

SINGLE_THREADED_SOLUTIONS = [
    'braid_multiplication_sequential_non_optimized',
    'braid_multiplication_sequential_memory',
    'braid_multiplication_sequential_precompute'
]
MULTI_THREADED_SOLUTIONS = [
    'braid_multiplication_parallel',
]

SOLUTIONS_FOLDER = 'braid_mults_solutions'  # where we put our ready to run implementations

CSV_FILE = 'braid_mul.csv'

REPEATS = 10
default_logger = Logger('log.txt')
build_logger = Logger('build_log.txt')

if __name__ == '__main__':
    arg_parser = ArgumentParser()
    arg_parser.add_argument('depth', type=int, help='Depth for level paralellism')
    args = arg_parser.parse_args()
    runners = compile_programs(SINGLE_THREADED_SOLUTIONS, MULTI_THREADED_SOLUTIONS, 0, args.depth + 1, BraidMulRunner,
                               build_logger, SOLUTIONS_FOLDER, CXX_COMPILER_PATH)

    n = [1000, 5000, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000,150000,200000, 250000,
         300000,350000,400000,450000,500000,550000,600000,650000,700000,750000,800000,850000,900000,950000,1000000,
         2000000,  5000000, 7500000, 10000000]
    tests = [BraidMulTest(str(x), '42') for x in n]

    run_tests(runners, tests, REPEATS, CSV_FILE, default_logger, False)
