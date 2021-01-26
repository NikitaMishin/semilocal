import subprocess
from utils.results import Result, BraidResult, CombingResult, LCSResult


class Runner:
    TIMEOUT = 60 * 20  ## in seconds

    def __init__(self, path, **kwargs):
        self.path = path
        self.name = path.split('/')[-1]

    def run(self, **kwargs) -> Result:
        pass


class BraidMulRunner(Runner):
    def __init__(self, depth, path, **kwargs):
        super().__init__(path, **kwargs)
        self.name += f'{depth}-depth'
        self.depth = depth

    def run(self, matrix_size, seed, **kwargs):
        try:
            compl_proc = subprocess.run([self.path, self.depth, matrix_size, seed],
                                        timeout=self.TIMEOUT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if compl_proc.stderr != b'':
                raise Exception(compl_proc.stderr)

            preprocess, algo, hash, n, seed = compl_proc.stdout.decode().split('\n')
            return BraidResult(int(preprocess), int(algo), int(hash), int(matrix_size), int(seed))
        except subprocess.TimeoutExpired:
            raise Exception('timeout')


class CombingRunner(Runner):
    def __init__(self, num_thd, path, **kwargs):
        super().__init__(path, **kwargs)
        self.name += f'{num_thd}-threaded'
        self.num_thread = num_thd

    def run(self, seq_a, seq_b, **kwargs):
        try:
            compl_proc = subprocess.run([self.path, self.num_thread, seq_a, seq_b],
                                        timeout=self.TIMEOUT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if compl_proc.stderr != b'':
                raise Exception(compl_proc.stderr)

            preprocess, algo, hash, a_size, b_size, a_name, b_name = compl_proc.stdout.decode().split('\n')
            return CombingResult(int(preprocess), int(algo), int(hash), int(a_size), int(b_size), a_name, b_name)
        except subprocess.TimeoutExpired:
            raise Exception('timeout')


class LcsRunner(Runner):
    def __init__(self, num_thd, path, **kwargs):
        super().__init__(path, **kwargs)
        self.name += f'{num_thd}-threaded'
        self.num_thread = num_thd

    def run(self, seq_a, seq_b, **kwargs):
        try:
            compl_proc = subprocess.run([self.path, self.num_thread, seq_a, seq_b],
                                        timeout=self.TIMEOUT, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            if compl_proc.stderr != b'':
                raise Exception(compl_proc.stderr)

            preprocess, algo, score, a_size, b_size, a_name, b_name = compl_proc.stdout.decode().split('\n')
            return LCSResult(int(preprocess), int(algo), int(score), int(a_size), int(b_size), a_name, b_name)
        except subprocess.TimeoutExpired:
            raise Exception('timeout')
