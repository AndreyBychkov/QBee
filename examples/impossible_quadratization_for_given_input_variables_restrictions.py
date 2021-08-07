import sympy as sp
from qbee import *
from functools import partial
from time import time

if __name__ == '__main__':
    c1, c2, c3, c4 = sp.symbols("c1, c2, c3, c4")
    A, Ea, Ru = sp.symbols("A, Ea, Ru")
    T = sp.Symbol("T")
    dc1, dc2, dc3, dc4 = derivatives([c1, c2, c3, c4])
    tmp = -A * sp.exp(-Ea / (Ru * T)) * c1 ** 0.2 * c2 ** 1.3
    system = EquationSystem([
        sp.Eq(dc1, tmp),
        sp.Eq(dc2, 2 * tmp),
        sp.Eq(dc3, -tmp),
        sp.Eq(dc4, -2 * tmp)],
        parameter_variables=[A, Ea, Ru],
        input_variables=[T])

    print("Original nonlinear system:")
    system.print()
    print("=" * 50)
    timeout = partial(pruning_by_elapsed_time, start_t=time(), max_t=10)
    # {T: 1} means than T can have a derivative of order at most one => T'
    quadr_system = polynomialize_and_quadratize(system, {T: 1}, pruning_functions=[timeout] + default_pruning_rules)
    if quadr_system:
        print(quadr_system)