import sympy as sp
from qbee import *


if __name__ == '__main__':
    x1, x2 = functions("x1, x2")
    system = [
        (x1, x1**3 + x2**2),
        (x2, x1 + x2)
    ]
    quadr_system = quadratize(system)
    print(quadr_system)
