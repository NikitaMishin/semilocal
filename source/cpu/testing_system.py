import itertools
import os
import subprocess
from argparse import ArgumentParser
from os import listdir, mkdir
from os.path import join, exists
from statistics import mean, stdev
from subprocess import run as shell_run
from subprocess import PIPE
from typing import List, Tuple, Dict
import fastaparser
from random import choice


# from tqdm import tqdm
# import numpy as np


class Logger:
    def __init__(self, log_file):
        self.file = open(log_file, 'wt', 1)

    def write(self, msg):
        print(msg)
        self.file.write(msg + '\n')

    def close(self):
        self.file.close()


LOG_FILE = Logger('log.txt')
BUILD_LOG_FILE = Logger('build_log.txt')
SOLUTIONS_FOLDER = 'semilocal_solutions'
TEST_FOLDER = 'semilocal_tests'
REAL_FASTA = 'semilocal_real_fasta'
SYNTHETIC_FASTA = 'synthetic_fasta'
ALPHABET_SIZE = [2, 4, 26, 33]

RESULT_CSV_FILE = 'results.csv'

compile_commands = [
                       ['rm', '-rf', SOLUTIONS_FOLDER],
                       ['rm', '-rf', TEST_FOLDER],
                       ['mkdir', SOLUTIONS_FOLDER],
                       ['mkdir', TEST_FOLDER],
                       ['mkdir', os.path.join(TEST_FOLDER, REAL_FASTA)],
                       ['mkdir', os.path.join(TEST_FOLDER, SYNTHETIC_FASTA)],
                       ['cmake', f'-B{SOLUTIONS_FOLDER}', '-H.', '-DCMAKE_CONFIGURATION_TYPES=Release'],
                       ['make', '-C', SOLUTIONS_FOLDER],
                   ] + [['mkdir', os.path.join(TEST_FOLDER, SYNTHETIC_FASTA, str(size))] for size in ALPHABET_SIZE]

REPEATS = 30


class Test:
    DELIM = '||||'

    def __init__(self, sequence_id_1: str, size_id_1: str, sequence_id_2: str, size_id_2: str):
        self.sequence_id_1 = sequence_id_1
        self.sequence_id_2 = sequence_id_2
        self.size_1 = size_id_1
        self.size_2 = size_id_2

    @property
    def name(self):
        return self.sequence_id_1 + self.DELIM + self.sequence_id_2


class Result:
    def __init__(self, elapsed_time_preprocess: int, elapsed_time_algo: int, score: int):
        self.score = score
        self.elapsed_time_algo = elapsed_time_algo
        self.elapsed_time_preprocess = elapsed_time_preprocess


class Runner:
    TIMEOUT = 60 * 30  ## in second

    def __init__(self, path, **kwargs):
        self.path = path
        self.name = path.split('/')[-1]

    def run(self, fasta_file_a, fasta_file_b) -> Result:
        pass


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
    def __init__(self, path, num_threads, **kwargs):
        super().__init__(path, **kwargs)
        self.name += f'{num_threads}-threaded'
        self.num_threads = num_threads

    def run(self, fasta_file_a, fasta_file_b) -> Result:
        try:
            compl_proc = subprocess.run([self.path, self.num_threads, fasta_file_a, fasta_file_b],
                                        timeout=self.TIMEOUT, stdout=PIPE, stderr=PIPE)
            if compl_proc.stderr != b'':
                raise Exception(compl_proc.stderr)
            preprocess, algo, score = compl_proc.stdout.decode().split('\n')
            return Result(int(preprocess), int(algo), int(score))
        except subprocess.TimeoutExpired:
            raise Exception('timeout')


def compile_programs(specified_solutions_single_threaded, specified_solutions_multi_threaded, max_thds: int) -> List[
    Runner]:
    print('Building programs...')

    for command in compile_commands:
        comp_proc = shell_run(command, stdout=BUILD_LOG_FILE.file, stderr=BUILD_LOG_FILE.file)
        if comp_proc.returncode != 0:
            exit(0)
    exec = [file for file in listdir(SOLUTIONS_FOLDER) if 'make' not in file.lower()]

    print(f'Built C++ solutions: {", ".join(exec)}')
    runners = \
        [CRunner(join(SOLUTIONS_FOLDER, ex), str(thds))
         for thds in range(1, max_thds + 1) for ex in exec if ex in specified_solutions_multi_threaded]

    runners += [CRunner(join(SOLUTIONS_FOLDER, ex), '1') for ex in exec if ex in specified_solutions_single_threaded]
    print(f'Runners  C++ to execute: {", ".join(exec)}')
    return runners


#
#
#
def build_tests(folder_with_fasta_files: str, only_generated=False) -> Tuple[Tuple[str, int], Tuple[str, int]]:
    with open("fasta_file.fasta") as fasta_file:
        parser = fastaparser.Reader(fasta_file)
        for seq in parser:
            # seq is a FastaSequence object
            print('ID:', seq.id)
            print('Description:', seq.description)
            print('Sequence:', seq.sequence_as_string())
            print()


def update_csv(runners: List[Runner], results: Dict[str, Dict[str, List[Result]]]):
    with open(RESULT_CSV_FILE, 'w') as f:
        header = 'Test name'
        for runner in runners:
            header += f',{runner.name}_time_preprocess_mean,{runner.name}_time_preprocess_std' + \
                      f',{runner.name}_time_algo_mean,{runner.name}_time_algo_std' + \
                      f',{runner.name}_score'
        f.write(header + '\n')
        for test, measures in results.items():
            line = test
            for runner in runners:
                if len(measures[runner.name]) > 0:
                    elapsed_time_preprocess = list(
                        map(lambda res: int(res.elapsed_time_preprocess), measures[runner.name]))
                    elapsed_time_algo = list(map(lambda res: int(res.elapsed_time_algo), measures[runner.name]))

                    mean_algo = mean(elapsed_time_algo)
                    mean_preprocessed = mean(elapsed_time_preprocess)
                    std_algo = stdev(elapsed_time_algo) if len(elapsed_time_algo) >= 2 else '-'
                    std_preprocessed = stdev(elapsed_time_preprocess) if len(elapsed_time_preprocess) >= 2 else '-'
                    prefix_score = measures[runner.name][0].score
                else:
                    mean_algo = std_algo = mean_preprocessed = std_preprocessed = prefix_score = '-'
                line += f',{mean_preprocessed},{std_preprocessed},{mean_algo},{std_algo},{prefix_score}'
            f.write(line + '\n')


# a = '/home/nikita/projects/semilocal/source/cpu/bug.fasta'

def run_tests(runners: List[Runner], tests: List[Test], repeats: int):
    run_strategy = RunStrategy(runners, tests, repeats)
    results = {
        test.name: {
            runner.name: [] for runner in runners
        } for test in tests
    }

    update_csv(runners, results)

    info_log = {}

    for runner, test in run_strategy.strategy:

        if isinstance(runner, str) and runner == 'iteration':
            update_csv(runners, results)
            continue

        if info_log.get((runner.name, test.name), '') == 'failed':
            LOG_FILE.write(f'Skip test {test.name} for {runner.name}  due to previous failure')
            continue

        LOG_FILE.write(f'{runner.name} work on {test.name} test')
        try:
            result = runner.run(test.sequence_id_1, test.sequence_id_2)
            results[test.name][runner.name].append(result)
            LOG_FILE.write(
                f'{runner.name} complete {test.name} in {result.elapsed_time_algo + result.elapsed_time_preprocess} ms')
        except Exception as e:
            LOG_FILE.write(f'{runner.name} failed on test {test.name} because of {e}')
            info_log[(runner.name, test.name)] = 'failed'
    update_csv(runners, results)


def generate_sequence(destination_folder: str, name: str, size: int, alphabet_size: int) -> Tuple[str, int]:
    path = os.path.join(destination_folder, f'generated_{name}')
    with open(path, 'w+') as f:
        seq = ''.join([str(choice(range(alphabet_size))) for _ in range(size)])
        idx = f'generated_{name}'
        f.writelines([idx + '\n', seq + '\n'])
        return path, size


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
    LOG_FILE.write(f'Built real_test set\nTotal pairs is {len(tests)}')
    return list(map(lambda x: Test(*x[0], *x[1]), tests))


def build_synthetic_test_set(alphabet_sizes, seq_sizes: List[int]) -> List[Test]:
    cnt = 0
    syntethic_test_set = []
    size = 0
    for alphabet_size in alphabet_sizes:
        test_set = []
        for seq_size in seq_sizes * 2:
            test_set.append(
                generate_sequence(os.path.join(TEST_FOLDER, SYNTHETIC_FASTA, str(alphabet_size)), str(cnt), seq_size,
                                  alphabet_size))
            cnt = cnt + 1
        tests = list(itertools.combinations(test_set, 2))
        size += len(tests)
        syntethic_test_set.extend(tests)
    LOG_FILE.write(f'Built synthetic_test_set\nTotal pairs is {size}')
    return list(map(lambda x: Test(*x[0], *x[1]), syntethic_test_set))


# if(True):
# runners = compile_programs(todo)

# build_real_test_set(
#
# )
# compile
# build
# build
# run
class TestingSystem:
    def run(self):
        pass


# class CTestingSystem()

runners = compile_programs([], ['prefix_lcs'], 4)
tests = build_real_test_set('/') + build_synthetic_test_set([33], [1280, 2560, 5120])
run_tests(runners, tests, 5)
