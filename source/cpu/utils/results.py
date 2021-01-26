class Result:
    pass


class CombingResult(Result):
    def __init__(self, elapsed_time_preprocess: int, elapsed_time_algo: int, hash: int, size_a: int, size_b: int,
                 name_a, name_b):
        self.hash = hash
        self.elapsed_time_algo = elapsed_time_algo
        self.elapsed_time_preprocess = elapsed_time_preprocess
        self.size_a = size_a
        self.size_b = size_b
        self.name_a = name_a
        self.name_b = name_b


class BraidResult(Result):
    def __init__(self, elapsed_time_preprocess, elapsed_time_algo, hash, n, seed):
        self.hash = hash
        self.elapsed_time_algo = elapsed_time_algo
        self.elapsed_time_preprocess = elapsed_time_preprocess
        self.n = n
        self.seed = seed


class LCSResult(Result):
    def __init__(self, elapsed_time_preprocess: int, elapsed_time_algo: int, score: int, size_a: int, size_b: int,
                 name_a, name_b):
        self.score = score
        self.elapsed_time_algo = elapsed_time_algo
        self.elapsed_time_preprocess = elapsed_time_preprocess
        self.size_a = size_a
        self.size_b = size_b
        self.name_a = name_a
        self.name_b = name_b
