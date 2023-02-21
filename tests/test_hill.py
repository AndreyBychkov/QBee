import pytest
from qbee import *
from examples import generate_hill


strategies = [default_scoring, aeqd_scoring, smd_scoring]
pruning_funcs = [
    [pruning_by_quadratic_upper_bound],
    [pruning_by_squarefree_graphs],
    [pruning_by_squarefree_graphs, pruning_by_quadratic_upper_bound]
]
systems = dict([(d, generate_hill(d)) for d in [3, 6, 8, 10, 15]])
order = systems.keys()

@pytest.mark.benchmark
@pytest.mark.parametrize('strat', strategies)
@pytest.mark.parametrize('pruning', pruning_funcs)
@pytest.mark.parametrize('ord', order)
def test_hill(benchmark, ord, strat, pruning):
    benchmark(quadratize, systems[ord], (), strat, True, pruning)
