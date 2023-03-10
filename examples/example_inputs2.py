import sympy as sp
from qbee import *

if __name__ == '__main__':
    x, u = functions("x, u")
    system = [
        (x, x**2 * u),
    ]
    quadr_system = quadratize(system)
    quadr_system.print()
