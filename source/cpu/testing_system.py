import itertools
import os
import subprocess
from argparse import ArgumentParser
from os import listdir
from os.path import join
from random import choice, shuffle, seed, choices, uniform
from statistics import mean, stdev
from subprocess import PIPE
from subprocess import run as shell_run
from typing import List, Tuple, Dict

import fastaparser

CXX_COMPILER_PATH = '/usr/bin/g++-10'
SINGLE_THREADED_SOLUTIONS = ['prefix_lcs', 'prefix_lcs_via_braid_bits_4symbol_int',
                             'prefix_lcs_via_braid_bits_4symbol_long_long',
                             'prefix_lcs_via_braid_bits_binary_int',
                             'prefix_lcs_via_braid_bits_binary_long_long',
                             'semi_local_lcs_sequential',
                             'prefix_lcs_via_braid_sequential'
                             ]
MULTI_THREADED_SOLUTIONS = ['prefix_lcs_via_braid_bits_4symbol_mpi_int',
                            'prefix_lcs_via_braid_bits_4symbol_mpi_long_long',
                            'prefix_lcs_via_braid_bits_binary_mpi_int',
                            'prefix_lcs_via_braid_bits_binary_mpi_long_long',
                            'semi_local_lcs_mpi'
                            ]
MAX_THREADS = 4
REPEATS = 5
SOLUTIONS_FOLDER = 'semilocal_solutions'
TEST_FOLDER = 'semilocal_tests'
REAL_FASTA = 'semilocal_real_fasta'
SYNTHETIC_FASTA = 'synthetic_fasta'
ALPHABETS = ['ac', 'acgt']
SYNTHETHIC_SIZES = [1024, 1024, 1024 * 2, 1024 * 2,
                    # 1024 * 4, 1024 * 4, 1024 * 8, 1024 * 8, 1024 * 16, 1024 * 16,
                    # 1024 * 32, 1024 * 32,
                    # 1024 * 64, 1024 * 64, 1024 * 128, 1024 * 128
                    ]
NO_RESTRICT = list(map(len, ALPHABETS)) + [4]

solution_to_support_alphabet_resticted = {
    'prefix_lcs_via_braid_bits_4symbol_int': [4,2],
    'prefix_lcs_via_braid_bits_4symbol_long_long': [4,2],
    'prefix_lcs_via_braid_bits_4symbol_mpi_int': [4,2],
    'prefix_lcs_via_braid_bits_4symbol_mpi_long_long': [4,2],
    'prefix_lcs_via_braid_bits_binary_int': [2],
    'prefix_lcs_via_braid_bits_binary_long_long': [2],
    'prefix_lcs_via_braid_bits_binary_mpi_int': [2],
    'prefix_lcs_via_braid_bits_binary_mpi_long_long': [2],
}


class Logger:
    def __init__(self, log_file):
        self.file = open(log_file, 'wt', 1)

    def write(self, msg):
        # print(msg)
        self.file.write(msg + '\n')

    def close(self):
        self.file.close()


default_logger = Logger('log.txt')
build_logger = Logger('build_log.txt')


def get_compile_commands(flags: Dict[str, bool]):
    rm_and_mkdirs = [['rm', '-rf', SOLUTIONS_FOLDER],
                     ['rm', '-rf', TEST_FOLDER],
                     ['mkdir', SOLUTIONS_FOLDER],
                     ['mkdir', TEST_FOLDER],
                     ['mkdir', os.path.join(TEST_FOLDER, REAL_FASTA)],
                     ['mkdir', os.path.join(TEST_FOLDER, SYNTHETIC_FASTA)]] + [
                        ['mkdir', os.path.join(TEST_FOLDER, SYNTHETIC_FASTA, str(size))] for size in
                        map(len, ALPHABETS)]

    cmake_command = ['cmake', f'-B{SOLUTIONS_FOLDER}',
                     '-H.', '-DCMAKE_BUILD_TYPE=Release',
                     f'-DCMAKE_CXX_COMPILER={CXX_COMPILER_PATH}']
    if flags.get('-fno-tree-vectorize', False):
        cmake_command += '-DCMAKE_CXX_FLAGS= -fno-tree-vectorize'

    make = ['make', '-C', SOLUTIONS_FOLDER]
    rm_and_mkdirs.extend([cmake_command, make])

    return rm_and_mkdirs


class Test:
    DELIM = '||||'

    def __init__(self, sequence_id_1: str, size_id_1: str, sequence_id_2: str, size_id_2: str, alphabet_size: int):
        self.sequence_id_1 = sequence_id_1
        self.sequence_id_2 = sequence_id_2
        self.size_1 = size_id_1
        self.size_2 = size_id_2
        self.alphabet_size = alphabet_size

    @property
    def name(self):
        return self.sequence_id_1 + self.DELIM + self.sequence_id_2


class Result:
    def __init__(self, elapsed_time_preprocess: int, elapsed_time_algo: int, score: int, size_a: int, size_b: int):
        self.score = score
        self.elapsed_time_algo = elapsed_time_algo
        self.elapsed_time_preprocess = elapsed_time_preprocess
        self.size_a = size_a
        self.size_b = size_b


class Runner:
    TIMEOUT = 60 * 30  ## in second

    def __init__(self, path, **kwargs):
        self.path = path
        self.name = path.split('/')[-1]

    def run(self, fasta_file_a, fasta_file_b) -> Result:
        pass

    def supported_alphabet(self) -> List[int]:
        pass
        # return self.supported_alphabet_sizes


class RunStrategy:

    def __init__(self, runners: List[Runner], tests: List[Test], repeats, method='simple'):
        self.runners = runners
        self.tests = tests
        self.method = method
        self.repeats = repeats
        if self.method == 'simple':
            self.strategy = self._simple_strategy()
        else:
            print('Provide correct strategy for build order of testing')

    def _simple_strategy(self) -> List[Tuple[Runner, Test]]:
        return [
            (runner, test)
            for _ in range(self.repeats)
            for test in self.tests
            for runner in self.runners + ['iteration']
        ]


class CRunner(Runner):
    def __init__(self, path, num_threads, supported_alphabet_sizes, **kwargs):
        super().__init__(path, **kwargs)
        self.name += f'{num_threads}-threaded'
        self.num_threads = num_threads
        self.supported_alphabet_sizes = supported_alphabet_sizes

    def run(self, fasta_file_a, fasta_file_b) -> Result:
        try:
            compl_proc = subprocess.run([self.path, self.num_threads, fasta_file_a, fasta_file_b],
                                        timeout=self.TIMEOUT, stdout=PIPE, stderr=PIPE)
            if compl_proc.stderr != b'':
                raise Exception(compl_proc.stderr)
            preprocess, algo, score, size_a, size_b = compl_proc.stdout.decode().split('\n')
            return Result(int(preprocess), int(algo), int(score), int(size_a), int(size_b))
        except subprocess.TimeoutExpired:
            raise Exception('timeout')

    def supported_alphabet(self) -> List[int]:
        return self.supported_alphabet_sizes


def compile_programs(specified_solutions_single_threaded, specified_solutions_multi_threaded, max_thds: int,
                     flags: Dict[str, bool]) -> List[Runner]:
    print('Building programs...')

    for command in get_compile_commands(flags):
        comp_proc = shell_run(command, stdout=build_logger.file, stderr=build_logger.file)
        if comp_proc.returncode != 0:
            exit(0)
    exec = [file for file in listdir(SOLUTIONS_FOLDER) if 'make' not in file.lower()]

    print(f'Built C++ solutions: {", ".join(exec)}')
    runners = \
        [CRunner(join(SOLUTIONS_FOLDER, ex), str(thds), solution_to_support_alphabet_resticted.get(ex, NO_RESTRICT))
         for thds in range(1, max_thds + 1) for ex in exec if ex in specified_solutions_multi_threaded]

    runners += [CRunner(join(SOLUTIONS_FOLDER, ex), '1', solution_to_support_alphabet_resticted.get(ex, NO_RESTRICT))
                for ex in exec if ex in specified_solutions_single_threaded]
    print(f'Runners  C++ to execute: {", ".join(map(lambda x: x.name, runners))}')
    return runners


def update_csv(runners: List[Runner], results: Dict[str, Dict[str, List[Result]]], flags):
    with open(flags['result_file'], 'w') as f:
        header = 'Test name, seq_a_size, seq_b_size'
        for runner in runners:
            header += f',{runner.name}_time_preprocess_mean,{runner.name}_time_preprocess_std' + \
                      f',{runner.name}_time_algo_mean,{runner.name}_time_algo_std' + \
                      f',{runner.name}_score'
        f.write(header + '\n')
        for test, measures in results.items():
            line = ''
            seq_a_size = 0
            seq_b_size = 0

            for runner in runners:
                if len(measures[runner.name]) > 0:
                    elapsed_time_preprocess = list(
                        map(lambda res: int(res.elapsed_time_preprocess), measures[runner.name]))
                    elapsed_time_algo = list(map(lambda res: int(res.elapsed_time_algo), measures[runner.name]))
                    seq_a_size = measures[runner.name][0].size_a
                    seq_b_size = measures[runner.name][0].size_b
                    mean_algo = mean(elapsed_time_algo)
                    mean_preprocessed = mean(elapsed_time_preprocess)
                    std_algo = stdev(elapsed_time_algo) if len(elapsed_time_algo) >= 2 else '-'
                    std_preprocessed = stdev(elapsed_time_preprocess) if len(elapsed_time_preprocess) >= 2 else '-'
                    prefix_score = measures[runner.name][0].score
                else:
                    mean_algo = std_algo = mean_preprocessed = std_preprocessed = prefix_score = seq_a_size = seq_b_size = '-'
                line += f',{mean_preprocessed},{std_preprocessed},{mean_algo},{std_algo},{prefix_score}'
            f.write( test + f',{seq_a_size},{seq_b_size}'+line + '\n')


# a = '/home/nikita/projects/semilocal/source/cpu/bug.fasta'

def run_tests(runners: List[Runner], tests: List[Test], repeats: int, flags):
    run_strategy = RunStrategy(runners, tests, repeats)
    results = {
        test.name: {
            runner.name: [] for runner in runners
        } for test in tests
    }

    update_csv(runners, results, flags)

    info_log = {}
    cnt = 0
    for runner, test in run_strategy.strategy:
        print(f'------{cnt} out of {len(run_strategy.strategy)}  ------')
        cnt += 1

        if isinstance(runner, str) and runner == 'iteration':
            update_csv(runners, results, flags)
            continue

        if info_log.get((runner.name, test.name), '') == 'failed':
            default_logger.write(f'Skip test {test.name} for {runner.name}  due to previous failure')
            continue
        if info_log.get((runner.name, test.name), '') == 'not support':
            # already logged that not support alphabet
            continue

        default_logger.write(f'{runner.name} work on {test.name} test')
        try:
            if test.alphabet_size in runner.supported_alphabet():
                result = runner.run(test.sequence_id_1, test.sequence_id_2)
                results[test.name][runner.name].append(result)
                default_logger.write(
                    f'{runner.name} complete {test.name} in {result.elapsed_time_algo + result.elapsed_time_preprocess} ms')
            else:
                default_logger.write(
                    f'{runner.name} not supported  test {test.name} with alphabet {test.alphabet_size}')
                info_log[(runner.name, test.name)] = 'not support'
        except Exception as e:
            default_logger.write(f'{runner.name} failed on test {test.name} because of {e}')
            info_log[(runner.name, test.name)] = 'failed'

    update_csv(runners, results, flags)


def generate_sequence(destination_folder: str, name: str, size: int, alphabet) -> Tuple[str, int, int]:
    seed(a=name)
    path = os.path.join(destination_folder, f'generated_{name}')
    weights = [uniform(0, 10) for _ in alphabet]
    with open(path, 'w+') as f:
        # lst = [str(choice(range(alphabet_size))) for _ in range(size)]
        lst = choices(list(alphabet), weights=weights, k=size)
        shuffle(lst)
        seq = ''.join(lst)
        idx = f'generated_{name}'
        f.writelines([idx + '\n', seq + '\n'])
        return path, size, len(alphabet)


def fasta_format_to_simple_format(fasta_file_path: str, new_dir_path) -> List[Tuple[str, int]]:
    """
    Convert fasta format with sequences to simple format of kind: name\n
                                                                  numerical_sequences\n
    One file per sequence
    :param fasta_file_path:
    :param new_dir_path:
    :return:
    """
    mapper = {
        'a': 0, 'c': 1, 'g': 2, 't': 3
    }
    files = []
    with open(fasta_file_path) as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            name = os.path.join(new_dir_path, seq.id)
            with open(name, 'w+') as new_simple_format:
                int_seq = ''.join(
                    list(map(lambda x: str(mapper.get(x, choice('0123'))), seq.sequence_as_string().lower())))
                new_simple_format.writelines([seq.id + '\n', int_seq + '\n'])
                files.append((name, len(int_seq)))
    return files


def build_real_test_set(path_to_fasta_files) -> List[Test]:
    fasta_files = []
    for root, dirs, files in os.walk(path_to_fasta_files):
        for file in files:
            if file.endswith('.fasta'):
                print(os.path.join(root, file))
                fasta_files.append(os.path.join(root, file))
    fasta_files = \
        list(map(lambda file: fasta_format_to_simple_format(file, os.path.join(TEST_FOLDER, REAL_FASTA)), fasta_files))
    path_size_tuples = []
    for lst in fasta_files:
        for path_size in lst:
            path_size_tuples.append(path_size)
    tests = list(itertools.combinations(path_size_tuples, 2))
    default_logger.write(f'Built real_test set\nTotal pairs is {len(tests)}')
    return list(map(lambda x: Test(x[0][0], x[0][1], x[1][0], x[1][1], 4), tests))


def build_synthetic_test_set(alphabets, seq_sizes: List[int]) -> List[Test]:
    cnt = 0
    syntethic_test_set = []
    size = 0
    for alphabet in alphabets:
        test_set = []
        for seq_size in seq_sizes:
            test_set.append(
                generate_sequence(os.path.join(TEST_FOLDER, SYNTHETIC_FASTA, str(len(alphabet))), str(cnt), seq_size,
                                  alphabet))
            cnt = cnt + 1
        tests = list(itertools.combinations(test_set, 2))
        size += len(tests)
        syntethic_test_set.extend(tests)
    default_logger.write(f'Built synthetic_test_set\nTotal pairs is {size}')
    return list(map(lambda x: Test(x[0][0], x[0][1], *x[1]), syntethic_test_set))


# class CTestingSystem()
if __name__ == '__main__':
    flags = {}
    arg_parser = ArgumentParser()
    arg_parser.add_argument('tests', type=str, help='Path to directory with real fasta files')
    arg_parser.add_argument('--result_file', type=str, nargs='?', help='Path to file where csv results will be stored',
                            default='result.csv')
    arg_parser.add_argument('--vectorization', nargs='?', type=bool, default=False,
                            help='flag that set vectorization option for g++')
    args = arg_parser.parse_args()

    flags['vectorization'] = args.vectorization
    flags['result_file'] = args.result_file

    runners = compile_programs(SINGLE_THREADED_SOLUTIONS, MULTI_THREADED_SOLUTIONS, MAX_THREADS, flags)
    tests = build_real_test_set(args.tests)
    tests += build_synthetic_test_set(ALPHABETS, SYNTHETHIC_SIZES)
    run_tests(runners, tests, REPEATS, flags)
    build_logger.close()
    default_logger.close()
