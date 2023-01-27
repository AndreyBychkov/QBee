import pytest
from qbee import *
from examples import generate_x_plus_1

generation = [default_generation]
scoring = [default_scoring, aeqd_scoring, smd_scoring]
pruning_funcs = [
    [pruning_by_quadratic_upper_bound],
    [pruning_by_squarefree_graphs],
    [pruning_by_squarefree_graphs, pruning_by_quadratic_upper_bound]
]
systems = dict([(d, generate_x_plus_1(d)) for d in [6, 8, 10, 12]])
order = systems.keys()

@pytest.mark.benchmark
@pytest.mark.parametrize('gen', generation)
@pytest.mark.parametrize('scoring', scoring)
@pytest.mark.parametrize('pruning', pruning_funcs)
@pytest.mark.parametrize('ord', order)
def test_x_plus_1(benchmark, ord, gen, scoring, pruning):
    benchmark(quadratize, systems[ord], (), True, gen, scoring, pruning)
