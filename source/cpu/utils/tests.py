from typing import List, Tuple

from utils.runners import Runner

class Test:
    DELIM = '||||'

class CombingTest(Test):


    def __init__(self, sequence_id_1: str, sequence_id_2: str):
        self.sequence_id_1 = sequence_id_1
        self.sequence_id_2 = sequence_id_2

    @property
    def name(self):
        return self.sequence_id_1 + self.DELIM + self.sequence_id_2


class BraidMulTest(Test):

    def __init__(self, n, seed):
        self.n = n
        self.seed = seed

    @property
    def name(self):
        return self.n + self.DELIM + self.seed


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
