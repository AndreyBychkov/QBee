import pytest
import sympy as sp
from qbee import *
from qbee.heuristics import *

R, x = sp.ring(["x"], sp.QQ)

heuristics = [default_score, aeqd_score, smd_score]
systems = dict([(d, PolynomialSystem([(x + 1) ** d])) for d in range(6, 15)])
order = systems.keys()


def quadratize(heuristics: Heuristics, order):
    algo = BranchAndBound(systems[order], heuristics=heuristics)
    algo.quadratize()


@pytest.mark.parametrize('heur', heuristics)
@pytest.mark.parametrize('ord', order)
def test_sigmoid(benchmark, heur, ord):
    benchmark(quadratize, heuristics=heur, order=ord)

