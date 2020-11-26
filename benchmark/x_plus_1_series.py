import pytest
import sympy as sp
from qbee import *
from qbee.heuristics import *
from qbee.quadratization import termination_by_square_bound, termination_by_best_nvars, termination_by_C4_bound

R, x = sp.ring(["x"], sp.QQ)

heuristics = [default_score, aeqd_score, smd_score]
terminations = [termination_by_square_bound, termination_by_best_nvars, termination_by_C4_bound]
systems = dict([(d, PolynomialSystem([(x + 1) ** d])) for d in range(6, 9)])
order = systems.keys()


def quadratize(order, heuristics: Heuristics, termination):
    algo = BranchAndBound(systems[order], heuristics=heuristics, early_termination=[termination])
    algo.quadratize()

@pytest.mark.parametrize('term', terminations)
@pytest.mark.parametrize('heur', heuristics)
@pytest.mark.parametrize('ord', order)
def test_sigmoid(benchmark, ord, heur, term):
    benchmark(quadratize, order=ord, heuristics=heur, termination=term)

