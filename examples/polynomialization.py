import sympy as sp
from qbee import *


if __name__ == '__main__':
    x = functions("x")
    poly_system = polynomialize([(x, sp.exp(x) + sp.exp(2 * x))])
    if poly_system:
        poly_system.print()
