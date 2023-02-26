import sympy as sp
from qbee import *

if __name__ == '__main__':
    x, u = functions("x1, u")
    system = [
        (x, x**2 * u),
    ]
    quadr_system = quadratize(system)
    print(quadr_system)
