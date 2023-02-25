import sympy as sp
from qbee import *

x, u = functions("x1, u")
system = [
    (x, x**2 * u),
]
quadr_system = quadratize(system)
print(quadr_system)
