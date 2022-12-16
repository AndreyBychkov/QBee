import sympy as sp
from qbee import *


if __name__ == '__main__':
    psi, theta, Apsi, Atheta, u = functions("pse, theta, Apsi, Atheta, u")
    bpsi, btheta, b, c0, c1, c2, c3, B, BD = parameters("bpsi, btheta, b, c0, c1, c2, c3, B, BD")
    system = [
        (psi, Apsi + bpsi - B * psi * (c0 + c1 * theta + c2 * theta**2 + c3 * theta**3)),
        (theta, Atheta + btheta + b * u + BD * psi * (c0 + c1 * theta + c2 * theta**2 + c3 * theta**3))
    ]

    quadr_system = polynomialize_and_quadratize(system, input_der_orders={Apsi: 0, Atheta: 0, u: 0})
    if quadr_system:
        print(quadr_system)
