import pytest
from qbee import *
from examples import generate_hard
from functools import partial


scorings = [default_scoring, aeqd_scoring, smd_scoring]
pruning_funcs = [
    [pruning_by_quadratic_upper_bound],
    [pruning_by_squarefree_graphs],
    [pruning_by_squarefree_graphs, pruning_by_quadratic_upper_bound]
]
systems = dict([(d, generate_hard(d)) for d in [3, 4]])
order = systems.keys()

@pytest.mark.benchmark
@pytest.mark.parametrize('scoring', scorings)
@pytest.mark.parametrize('pruning', pruning_funcs)
@pytest.mark.parametrize('ord', order)
def test_hard(benchmark, ord, scoring, pruning):
    benchmark(partial(quadratize, systems[ord], scoring=scoring, pruning_functions=pruning))

