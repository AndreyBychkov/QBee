import sympy as sp
from qbee import *

if __name__ == '__main__':
    x1, x2, u = functions("x1, x2, u")
    system = [
        (x1, x1 + x1 * u),
        (x2, x1**2 * u)
    ]
    quadr_system = quadratize(system, input_free=True)
    print(quadr_system)
