import pytest
from functools import partial
from qbee import *


def test_dym_equations():
    x = functions("x")
    u = multivariable_functions("u", [x])
    rhs = u ** 3 * u.diff(x, 3)
    system = [
        (x, x),
        (u, rhs),
        (u.diff(x), rhs.diff(x)),
        (u.diff(x, 2), rhs.diff(x, 2))
    ]
    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize(system, pruning_functions=[*default_pruning_rules,
                                                                  no_quad_pruning
                                                                  ])
    assert res is not None
    print(res)
