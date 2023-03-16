from sympy import *
from qbee import *

if __name__ == '__main__':
    psi, theta, w1, Dpsi, Dtheta = functions("psi theta w1 Dpsi Dtheta")
    u = functions("u")
    B, D, gamma, bpsi, btheta, b = parameters("B D gamma bpsi btheta b")
    theta_diff = Dtheta + btheta + b * u + B * D * psi * w1
    system = [
        (psi, Dpsi + bpsi - D * psi * w1),
        (theta, theta_diff),
        (w1, gamma * w1 / theta**2 * theta_diff)
    ]
    new_vars = quadratize_dimension_agnostic(system)
    print(new_vars)
