from sympy import *
from qbee import *

if __name__ == '__main__':
    psi, theta, w1, Dpsi, Dtheta = functions("psi theta w1 Dpsi Dtheta")
    u = functions("u")
    B, D, gamma, bpsi, btheta = parameters("B   D gamma bpsi btheta")
    theta_diff = Dtheta + btheta * u + B * D * psi * w1
    system = [
        (psi, Dpsi + bpsi * u - D * psi * w1),
        (theta, theta_diff),
        (w1, gamma * w1 / theta**2 * theta_diff)
    ]
    new_vars = quadratize_dimension_agnostic(system, non_duplicated_vars=[u])
    print(new_vars)
