import sympy as sp
from qbee import *


if __name__ == '__main__':
    x1, x2, x3, x4 = functions("x1, x2, x3, x4")
    u = functions("u")
    A, E, R = parameters("A, E, R")
    x1_diff = -A * sp.exp(-E / (R * u)) * x1**0.2 * x2**1.3
    system = [
        (x1, x1_diff),
        (x2, 2 * x1_diff),
        (x3, -x1_diff),
        (x4, -2 * x1_diff)
    ]
    laurent_system = polynomialize(system)
    laurent_system.print()
    quadr_system = polynomialize_and_quadratize(system, input_der_orders={u: 2})
    quadr_system.print()
