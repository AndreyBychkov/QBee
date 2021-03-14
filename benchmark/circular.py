import pytest
import sympy as sp
from qbee import *
from qbee.selection import *
from qbee.quadratization import pruning_by_quadratic_upper_bound, pruning_by_best_nvars, pruning_by_squarefree_graphs
from qbee.examples import generate_circular

R, x = sp.ring(["x"], sp.QQ)

# strategies = [default_score, aeqd_score, smd_score]  # Does not really change performance by now
strategies = [aeqd_strategy]
pruning_funcs = [pruning_by_quadratic_upper_bound, pruning_by_best_nvars, pruning_by_squarefree_graphs]
systems = dict([(d, generate_circular(d)) for d in range(3, 6)])
order = systems.keys()


def quadratize(order, strategy: SelectionStrategy, pruning):
    algo = BranchAndBound(systems[order], strategy=strategy, pruning_funcs=[pruning])
    algo.quadratize()


@pytest.mark.parametrize('pruning', pruning_funcs)
@pytest.mark.parametrize('strat', strategies)
@pytest.mark.parametrize('ord', order)
def test_circular(benchmark, ord, strat, pruning):
    benchmark(quadratize, order=ord, strategy=strat, termination=pruning)
