from os import listdir
from os.path import join
from statistics import mean, stdev
from subprocess import run as shell_run

from utils.tests import RunStrategy


def update_csv(runners, results, csv_file, is_combing):
    with open(csv_file, 'w') as f:
        if is_combing:
            header = 'Test name, a_name, b_name, size_a, size_b'
        else:
            header = 'Test name, n, seed'


        for runner in runners:
            header += f',{runner.name}_time_preprocess_mean,{runner.name}_time_preprocess_std' + \
                      f',{runner.name}_time_algo_mean,{runner.name}_time_algo_std' + \
                      f',{runner.name}_hash'
        f.write(header + '\n')
        for test, measures in results.items():
            line = ''
            n = seed = '-'
            a_name = b_name = size_a = size_b = '-'
            for runner in runners:
                if len(measures[runner.name]) > 0:
                    elapsed_time_preprocess = list(
                        map(lambda res: int(res.elapsed_time_preprocess), measures[runner.name]))
                    elapsed_time_algo = list(map(lambda res: int(res.elapsed_time_algo), measures[runner.name]))

                    mean_algo = mean(elapsed_time_algo)
                    mean_preprocessed = mean(elapsed_time_preprocess)

                    if is_combing:
                        a_name = measures[runner.name][0].name_a
                        b_name = measures[runner.name][0].name_b
                        size_a, size_b = measures[runner.name][0].size_a, measures[runner.name][0].size_b
                    else:
                        n = measures[runner.name][0].n
                        seed = measures[runner.name][0].seed

                    std_algo = stdev(elapsed_time_algo) if len(elapsed_time_algo) >= 2 else '-'
                    std_preprocessed = stdev(elapsed_time_preprocess) if len(elapsed_time_preprocess) >= 2 else '-'
                    hash = measures[runner.name][0].hash
                else:
                    mean_algo = std_algo = mean_preprocessed = std_preprocessed = hash = '-'
                line += f',{mean_preprocessed},{std_preprocessed},{mean_algo},{std_algo},{hash}'
            if is_combing:
                f.write(test + f',{a_name},{b_name},{size_a},{size_b}' + line + '\n')
            else:
                f.write(test + f',{n},{seed}' + line + '\n')


def run_tests(runners, tests, repeats: int, csv_file, default_logger, is_combing = False):
    run_strategy = RunStrategy(runners, tests, repeats)
    results = {
        test.name: {
            runner.name: [] for runner in runners
        } for test in tests
    }

    update_csv(runners, results, csv_file,is_combing)

    info_log = {}
    cnt = 0
    for runner, test in run_strategy.strategy:
        print(f'------{cnt} out of {len(run_strategy.strategy)}  ------')
        cnt += 1

        if isinstance(runner, str) and runner == 'iteration':
            update_csv(runners, results, csv_file,is_combing)
            continue

        if info_log.get((runner.name, test.name), '') == 'failed':
            default_logger.write(f'Skip test {test.name} for {runner.name}  due to previous failure')
            continue

        default_logger.write(f'{runner.name} work on {test.name} test')

        try:

            if not is_combing:
                result = runner.run(test.n, test.seed)
            else:
                result = runner.run(test.sequence_id_1, test.sequence_id_2)

            results[test.name][runner.name].append(result)

            default_logger.write(
                f'{runner.name} complete {test.name} in {result.elapsed_time_algo + result.elapsed_time_preprocess} ms')

        except Exception as e:
            default_logger.write(f'{runner.name} failed on test {test.name} because of {e}')
            info_log[(runner.name, test.name)] = 'failed'

    update_csv(runners, results, csv_file,is_combing)


def build_algos(sequential_algos, parallel_algos, min_bound, max_bound, folder_with_impls, RunnerClass):
    """
    Build BraudMul runners
    :param sequential_algos:
    :param parallel_algos:
    :param max_depth:
    :param folder_with_impls:
    :return:
    """
    implementations = [file for file in listdir(folder_with_impls) if file.lower() in sequential_algos + parallel_algos]

    runners = [RunnerClass(str(min_bound), join(folder_with_impls, ex)) for ex in implementations if
               ex in sequential_algos + parallel_algos]
    runners += [RunnerClass(str(depth), join(folder_with_impls, ex)) for depth in range(min_bound + 1, max_bound) for ex in
                implementations if ex in parallel_algos]

    print(f'Runners  C++ to execute: {", ".join(map(lambda x: x.name, runners))}')
    return runners


def compile_programs(specified_solutions_single_threaded, specified_solutions_multi_threaded, min_bound, max_bound,
                     RunnerClass, build_logger, sol_folder, compiler_path):
    from utils.compile import construct_compile_commands

    print('Building programs...')

    for command in construct_compile_commands(sol_folder, compiler_path):
        comp_proc = shell_run(command, stdout=build_logger.file, stderr=build_logger.file)
        if comp_proc.returncode != 0:
            exit(0)

    return build_algos(specified_solutions_single_threaded, specified_solutions_multi_threaded, min_bound, max_bound,
                       sol_folder, RunnerClass)
