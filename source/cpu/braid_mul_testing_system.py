from argparse import ArgumentParser
from os import listdir
from os.path import join
from statistics import mean, stdev
from subprocess import run as shell_run
from typing import List, Dict

from utils.logger import Logger
from utils.results import BraidResult
from utils.runners import BraidMulRunner, Runner
from utils.tests import RunStrategy, BraidMulTest

# /usr/bin/g++-10
CXX_COMPILER_PATH = '/usr/bin/g++'

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

REPEATS = 30
default_logger = Logger('log.txt')
build_logger = Logger('build_log.txt')


def build_braid_mul_algorithms(sequential_algos, parallel_algos, max_depth, folder_with_impls):
    """
    Build BraudMul runners
    :param sequential_algos:
    :param parallel_algos:
    :param max_depth:
    :param folder_with_impls:
    :return:
    """
    implementations = [file for file in listdir(folder_with_impls) if file.lower() in sequential_algos + parallel_algos]

    runners = [BraidMulRunner('0', join(folder_with_impls, ex)) for ex in implementations if
               ex in sequential_algos + parallel_algos]
    runners += [BraidMulRunner(str(depth), join(folder_with_impls, ex)) for depth in range(1, max_depth) for ex in
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

    return build_braid_mul_algorithms(specified_solutions_single_threaded, specified_solutions_multi_threaded,
                                      max_depth, SOLUTIONS_FOLDER)


def update_csv(runners: List[BraidMulRunner], results: Dict[str, Dict[str, List[BraidResult]]]):
    with open(CSV_FILE, 'w') as f:
        header = 'Test name, n, seed'
        for runner in runners:
            header += f',{runner.name}_time_preprocess_mean,{runner.name}_time_preprocess_std' + \
                      f',{runner.name}_time_algo_mean,{runner.name}_time_algo_std' + \
                      f',{runner.name}_hash'
        f.write(header + '\n')
        for test, measures in results.items():
            line = ''
            n = seed = '-'
            for runner in runners:
                if len(measures[runner.name]) > 0:
                    elapsed_time_preprocess = list(
                        map(lambda res: int(res.elapsed_time_preprocess), measures[runner.name]))
                    elapsed_time_algo = list(map(lambda res: int(res.elapsed_time_algo), measures[runner.name]))

                    mean_algo = mean(elapsed_time_algo)
                    mean_preprocessed = mean(elapsed_time_preprocess)

                    n = measures[runner.name][0].n
                    seed = measures[runner.name][0].seed

                    std_algo = stdev(elapsed_time_algo) if len(elapsed_time_algo) >= 2 else '-'
                    std_preprocessed = stdev(elapsed_time_preprocess) if len(elapsed_time_preprocess) >= 2 else '-'
                    hash = measures[runner.name][0].hash
                else:
                    mean_algo = std_algo = mean_preprocessed = std_preprocessed = hash = '-'
                line += f',{mean_preprocessed},{std_preprocessed},{mean_algo},{std_algo},{hash}'
            f.write(test + f',{n},{seed}' + line + '\n')


def run_tests(runners: List[BraidMulRunner], tests: List[BraidMulTest], repeats: int):
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

            result = runner.run(test.n, test.seed)
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
    arg_parser.add_argument('depth', type=int, help='Depth for level paralellism')
    args = arg_parser.parse_args()
    runners = compile_programs(SINGLE_THREADED_SOLUTIONS, MULTI_THREADED_SOLUTIONS, args.depth)

    n = [1000, 10000, 20000, 40000, 80000, 160000, 320000,500000, 640000,1000000,2000000,4000000,5000000,7500000,10000000]
    tests = [BraidMulTest(str(x), '42') for x in n]

    run_tests(runners, tests, REPEATS)
    # build_logger.close()
    # default_logger.close()
