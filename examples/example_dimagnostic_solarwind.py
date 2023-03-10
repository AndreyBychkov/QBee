from sympy import *
from qbee import *

if __name__ == '__main__':
    u, v, Du, Dv = functions("u v Du Dv")
    system = [
        (v, Du + Dv),
        (u, Du / v + Dv / v)
    ]
    new_vars = quadratize_dimension_agnostic(system)
    print(new_vars)
