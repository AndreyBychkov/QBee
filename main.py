from qbee import *
from examples import *
from functools import partial

if __name__ == '__main__':
    x, a = functions("x, a")
    system = [
        (x, a * x ** 5),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize(system, pruning_functions=[no_quad_pruning, *default_pruning_rules],
                                       optimize_parameters={a})
    assert res is not None
    print(res)
