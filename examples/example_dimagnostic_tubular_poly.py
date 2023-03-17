from sympy import *
from qbee import *

if __name__ == '__main__':
    psi, theta, w1, Dpsi, Dtheta = functions("psi theta w1 Dpsi Dtheta")
    u = functions("u")
    B, D, gamma, bpsi, btheta, b = parameters("B D gamma bpsi btheta b")
    c0, c1, c2, c3 = parameters("c0 c1 c2 c3")
    poly = c0 + c1 * theta + c2 * theta**2 + c3 * theta**3
    system = [
        (psi, Dpsi + bpsi - D * psi * poly),
        (theta, Dtheta + btheta + b * u + B * D * psi * poly),
    ]
    new_vars = quadratize_dimension_agnostic(system)
