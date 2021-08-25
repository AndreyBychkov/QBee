import pytest
from qbee import *
from examples import generate_lifeware_conjecture


strategies = [default_strategy, aeqd_strategy, smd_strategy]
pruning_funcs = [
    [pruning_by_quadratic_upper_bound],
    [pruning_by_squarefree_graphs],
    [pruning_by_squarefree_graphs, pruning_by_quadratic_upper_bound]
]
systems = dict([(d, generate_lifeware_conjecture(d)) for d in [2, 3]])
order = systems.keys()


@pytest.mark.parametrize('strat', strategies)
@pytest.mark.parametrize('pruning', pruning_funcs)
@pytest.mark.parametrize('ord', order)
def test_lifeware_conjecture(benchmark, ord, strat, pruning):
    benchmark(quadratize, systems[ord], strat, pruning)