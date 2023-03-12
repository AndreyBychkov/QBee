from sympy import *
from qbee import *

if __name__ == '__main__':
    u, v, Du, Dv = functions("u v Du Dv")
    C1 = parameters("C1")
    system = [
        (v, Du + C1 * Dv),
        (u, Du / v + C1 * Dv / v)
    ]
    new_vars = quadratize_dimension_agnostic(system)
    print(new_vars)
