import pytest
from qbee import *
from examples import generate_lifeware_conjecture


strategies = [default_scoring, aeqd_scoring, smd_scoring]
pruning_funcs = [
    [pruning_by_quadratic_upper_bound],
    [pruning_by_squarefree_graphs],
    [pruning_by_squarefree_graphs, pruning_by_quadratic_upper_bound]
]
systems = dict([(d, generate_lifeware_conjecture(d)) for d in [2, 3]])
order = systems.keys()

@pytest.mark.benchmark
@pytest.mark.parametrize('strat', strategies)
@pytest.mark.parametrize('pruning', pruning_funcs)
@pytest.mark.parametrize('ord', order)
def test_lifeware_conjecture(benchmark, ord, strat, pruning):
    benchmark(quadratize, systems[ord], (), True, strat, pruning)
