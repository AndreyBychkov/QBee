from examples import *
from functools import partial

from qbee.experimental import polynomialize_and_quadratize_pde, multivariable_functions

if __name__ == '__main__':
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (x, x),
        (u, u.diff(x) ** 2),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize_pde(system, pruning_functions=[no_quad_pruning, *default_pruning_rules])
    assert res is not None
    print(res)

