import sympy as sp
from qbee import *


if __name__ == '__main__':
    v1, v2, v3, v1inv, v2inv, v3inv = functions("v1, v2, v3, v1inv, v2inv, v3inv")
    c = parameters("c")
    #vt = vx * vinv - c * vx
    v1t = (v1 + v2) * v1inv - c * (v1 + v2)
    v2t = (v1 + v2 + v3) * v2inv + c * (v1 + v2 + v3)
    v3t = (v2 + v3) * v3inv + c * (v2 + v3)
    system = [
        (v1, v1t),
        (v2, v2t),
        (v3, v3t),
        (v1inv, -v1t * v1inv**2),
        (v2inv, -v2t * v2inv**2),
        (v3inv, -v3t * v3inv**2)
    ]

    quadr_system = polynomialize_and_quadratize(system)
    if quadr_system:
        print(quadr_system)
