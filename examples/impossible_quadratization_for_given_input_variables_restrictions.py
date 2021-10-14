import sympy as sp
from qbee import *
from functools import partial
from time import time

if __name__ == '__main__':
    c1, c2, c3, c4, T = functions("c1, c2, c3, c4, T")
    A, Ea, Ru = parameters("A, Ea, Ru")
    eq1 = -A * sp.exp(-Ea / (Ru * T)) * c1 ** 0.2 * c2 ** 1.3
    system = [
        (c1, eq1),
        (c2, 2 * eq1),
        (c3, -eq1),
        (c4, -2 * eq1)
    ]

    timeout = partial(pruning_by_elapsed_time, start_t=time(), max_t=10)
    # {T: 1} means than T can have a derivative of order at most one => T'
    quadr_system = polynomialize_and_quadratize(system, input_der_orders={T: 1},
                                                pruning_functions=[timeout] + default_pruning_rules)
    if quadr_system:
        print(quadr_system)
