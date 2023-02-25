import sympy as sp
from qbee import *

x1, x2, u = functions("x1, x2, u")
system = [
    (x1, x1 + x1 * u),
    (x2, x1**2 * u)
]
quadr_system = quadratize(system, input_der_orders = {u: 0})
print(quadr_system)
