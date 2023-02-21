import sympy as sp
from qbee import *


if __name__ == '__main__':
    v1, v2, v3 = functions("v1, v2, v3")
    c = parameters("c")
    #vt = vx * vinv - c * vx
    v1t = (v1 + v2) * v1**(-1) - c * (v1 + v2)
    v2t = (v1 + v2 + v3) * v2**(-1) + c * (v1 + v2 + v3)
    v3t = (v2 + v3) * v3**(-1) + c * (v2 + v3)
    system = [
        (v1, v1t),
        (v2, v2t),
        (v3, v3t)
    ]

    quadr_system = polynomialize_and_quadratize(system)
    if quadr_system:
        print(quadr_system)
