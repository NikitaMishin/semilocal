import os
from argparse import ArgumentParser
from os import listdir
from os.path import join
from statistics import mean, stdev
from subprocess import run as shell_run
from typing import List, Dict

from utils.logger import Logger
from utils.results import CombingResult
from utils.runners import Runner, CombingRunner
from utils.tests import RunStrategy, CombingTest

CXX_COMPILER_PATH = '/usr/bin/g++-10'

SINGLE_THREADED_SOLUTIONS = [
    'prefix_lcs',

]

MULTI_THREADED_SOLUTIONS = [
    'semi_local_parallel_iterative',
    'semi_local_parallel_hybrid',
    'semi_local_parallel_1and3_combined_iterative',
    'semi_local_parallel_inversea_iterative',
    'semi_local_parallel_withoutif_iterative',
    'semi_local_parallel_iterative',
]

SOLUTIONS_FOLDER = 'combing_solutions'  # where we put our ready to run implementations

CSV_FILE = 'results.csv'

REPEATS = 10
default_logger = Logger('log.txt')
build_logger = Logger('build_log.txt')


def build_combing_algorithms(sequential_algos, parallel_algos, max_thds, folder_with_impls):
    """
    Build combingRunners
    :param sequential_algos:
    :param parallel_algos:
    :param max_thds:
    :param folder_with_impls:
    :return:
    """
    implementations = [file for file in listdir(folder_with_impls) if file.lower() in sequential_algos + parallel_algos]

    runners = [CombingRunner('1', join(folder_with_impls, ex)) for ex in implementations if
               ex in sequential_algos + parallel_algos]

    runners += [CombingRunner(str(thds), join(folder_with_impls, ex)) for thds in range(2, max_thds + 1) for ex in
                implementations if ex in parallel_algos]

    print(f'Runners  C++ to execute: {", ".join(map(lambda x: x.name, runners))}')
    return runners


def compile_programs(specified_solutions_single_threaded, specified_solutions_multi_threaded, max_depth) -> List[
    Runner]:
    from utils.compile import construct_compile_commands

    print('Building programs...')

    for command in construct_compile_commands(SOLUTIONS_FOLDER, CXX_COMPILER_PATH):
        comp_proc = shell_run(command, stdout=build_logger.file, stderr=build_logger.file)
        if comp_proc.returncode != 0:
            exit(0)

    return build_combing_algorithms(specified_solutions_single_threaded, specified_solutions_multi_threaded, max_depth,
                                    SOLUTIONS_FOLDER)


def update_csv(runners: List[CombingRunner], results: Dict[str, Dict[str, List[CombingResult]]]):
    with open(CSV_FILE, 'w') as f:
        header = 'Test name, a_name, b_name, size_a, size_b'
        for runner in runners:
            header += f',{runner.name}_time_preprocess_mean,{runner.name}_time_preprocess_std' + \
                      f',{runner.name}_time_algo_mean,{runner.name}_time_algo_std' + \
                      f',{runner.name}_hash'
        f.write(header + '\n')
        for test, measures in results.items():
            line = ''
            a_name = b_name = size_a = size_b = '-'
            for runner in runners:
                if len(measures[runner.name]) > 0:
                    elapsed_time_preprocess = list(
                        map(lambda res: int(res.elapsed_time_preprocess), measures[runner.name]))
                    elapsed_time_algo = list(map(lambda res: int(res.elapsed_time_algo), measures[runner.name]))

                    mean_algo = mean(elapsed_time_algo)
                    mean_preprocessed = mean(elapsed_time_preprocess)

                    name_a = measures[runner.name][0].name_a
                    name_b = measures[runner.name][0].name_b
                    size_a, size_b = measures[runner.name][0].size_a, measures[runner.name][0].size_b

                    std_algo = stdev(elapsed_time_algo) if len(elapsed_time_algo) >= 2 else '-'
                    std_preprocessed = stdev(elapsed_time_preprocess) if len(elapsed_time_preprocess) >= 2 else '-'
                    hash = measures[runner.name][0].hash
                else:
                    mean_algo = std_algo = mean_preprocessed = std_preprocessed = hash = '-'
                line += f',{mean_preprocessed},{std_preprocessed},{mean_algo},{std_algo},{hash}'
            f.write(test + f',{a_name},{b_name},{size_a},{size_b}' + line + '\n')


def run_tests(runners: List[CombingRunner], tests: List[CombingTest], repeats: int):
    run_strategy = RunStrategy(runners, tests, repeats)
    results = {
        test.name: {
            runner.name: [] for runner in runners
        } for test in tests
    }

    update_csv(runners, results)

    info_log = {}
    cnt = 0
    for runner, test in run_strategy.strategy:
        print(f'------{cnt} out of {len(run_strategy.strategy)}  ------')
        cnt += 1

        if isinstance(runner, str) and runner == 'iteration':
            update_csv(runners, results)
            continue

        if info_log.get((runner.name, test.name), '') == 'failed':
            default_logger.write(f'Skip test {test.name} for {runner.name}  due to previous failure')
            continue

        default_logger.write(f'{runner.name} work on {test.name} test')

        try:

            result = runner.run(test.sequence_id_1, test.sequence_id_2)
            results[test.name][runner.name].append(result)
            default_logger.write(
                f'{runner.name} complete {test.name} in {result.elapsed_time_algo + result.elapsed_time_preprocess} ms')

        except Exception as e:
            default_logger.write(f'{runner.name} failed on test {test.name} because of {e}')
            info_log[(runner.name, test.name)] = 'failed'

    update_csv(runners, results)


# class CTestingSystem()
if __name__ == '__main__':
    arg_parser = ArgumentParser()
    arg_parser.add_argument('max_thds', type=int, help='Threads')
    arg_parser.add_argument('tests', type=str, help='Path to dataset with *.fna files, our specific format')
    args = arg_parser.parse_args()

    runners = compile_programs(SINGLE_THREADED_SOLUTIONS, MULTI_THREADED_SOLUTIONS, args.max_thds)

    test_cases = []
    for name in listdir(args.tests):
        with open(os.path.join(args.tests,name), 'r') as f:
            _, size = f.readline(), int(f.readline())
            test_cases.append((size, os.path.join(args.tests, name)))

    tests = [CombingTest(x[1], y[1]) for x in test_cases for y in test_cases if x[1] != y[1] and x[0] == y[0] and x[0]]

    run_tests(runners, tests, REPEATS)
