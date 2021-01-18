import pytest
import sympy as sp
from qbee import *
from qbee.heuristics import *
from qbee.quadratization import termination_by_square_bound, termination_by_best_nvars, termination_by_C4_bound
from qbee.examples import generate_hill

R, x = sp.ring(["x"], sp.QQ)

heuristics = [default_score, aeqd_score, smd_score]
terminations = [termination_by_square_bound, termination_by_best_nvars, termination_by_C4_bound]
systems = dict([(d, generate_hill(d)) for d in [2, 3, 4, 5, 6, 7, 8, 10, 15]])
order = systems.keys()


def quadratize(order, heuristics: Heuristics, termination):
    algo = BranchAndBound(systems[order], heuristics=heuristics, early_termination=[termination])
    algo.quadratize()


@pytest.mark.parametrize('term', terminations)
@pytest.mark.parametrize('heur', heuristics)
@pytest.mark.parametrize('ord', order)
def test_hill(benchmark, ord, heur, term):
    benchmark(quadratize, order=ord, heuristics=heur, termination=term)