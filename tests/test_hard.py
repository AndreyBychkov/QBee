import pytest
from qbee import *
from examples import generate_hard

generation = [default_generation]
scoring = [default_scoring, aeqd_scoring, smd_scoring]
pruning_funcs = [
    [pruning_by_quadratic_upper_bound],
    [pruning_by_squarefree_graphs],
    [pruning_by_squarefree_graphs, pruning_by_quadratic_upper_bound]
]
systems = dict([(d, generate_hard(d)) for d in [3, 4]])
order = systems.keys()

@pytest.mark.benchmark
@pytest.mark.parametrize('gen', generation)
@pytest.mark.parametrize('scoring', scoring)
@pytest.mark.parametrize('pruning', pruning_funcs)
@pytest.mark.parametrize('ord', order)
def test_hard(benchmark, ord, gen, scoring, pruning):
    benchmark(quadratize, systems[ord], (), True, gen, scoring, pruning)
