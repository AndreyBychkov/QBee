import pytest
import sympy as sp
from qbee import *
from qbee.heuristics import *
from qbee.quadratization import pruning_by_quadratic_upper_bound, pruning_by_best_nvars, pruning_by_squarefree_graphs
from qbee.examples import generate_lifeware_conjecture

R, x = sp.ring(["x"], sp.QQ)

# heuristics = [default_score, aeqd_score, smd_score]  # Does not really change performance by now
heuristics = [default_score]
terminations = [pruning_by_quadratic_upper_bound, pruning_by_best_nvars, pruning_by_squarefree_graphs]
systems = dict([(d, generate_lifeware_conjecture(d)) for d in [2, 3]])
order = systems.keys()


def quadratize(order, heuristics: Heuristics, termination):
    algo = BranchAndBound(systems[order], heuristics=heuristics, early_termination=[termination])
    algo.quadratize()


@pytest.mark.parametrize('term', terminations)
@pytest.mark.parametrize('heur', heuristics)
@pytest.mark.parametrize('ord', order)
def test_long_monomials(benchmark, ord, heur, term):
    benchmark(quadratize, order=ord, heuristics=heur, termination=term)
