from sympy import *
from qbee import *

if __name__ == '__main__':
    x, Dx = functions("x Dx")
    system = [
        (x, x**2 * Dx)
    ]
    new_vars = quadratize_dimension_agnostic(system)
    print(new_vars)
