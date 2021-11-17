import sympy as sp
from qbee import *

def variables_repetitions_in_polynomilaization():
    n = 2
    psi = sp.symbols([f"psi{i}" for i in range(n)])
    theta = sp.symbols([f"theta{i}" for i in range(n)])

    bpsi = sp.symbols([f"bpsi{i}" for i in range(n)])
    btheta = sp.symbols([f"btheta{i}" for i in range(n)])
    Ap = sp.symbols([f"Ap_{i}_{j}" for i in range(n) for j in range(n)])
    At = sp.symbols([f"At_{i}_{j}" for i in range(n) for j in range(n)])
    B, BD, gamma = sp.symbols(["B", "BD", "gamma"])
    u = sp.symbols("u")

    dpsi = derivatives(psi)
    dtheta = derivatives(theta)

    tmp = [sp.exp(-gamma / theta[i]) for i in range(n)]

    eqs = []
    for i in range(n):
        eqs.append(
            sp.Eq(dpsi[i], sum([Ap[i * n + j] * psi[j] for j in range(n)]) + bpsi[i] * u - B * tmp[i] * psi[i])
        )
    for i in range(n):
        eqs.append(
            sp.Eq(dtheta[i], sum([At[i * n + j] * theta[j] for j in range(n)]) + btheta[i] * u + BD * tmp[i] * psi[i])
        )

    system = EquationSystem(
        eqs,
        parameter_variables=bpsi + btheta + [B, BD, gamma] + Ap + At,
        input_variables=[u])
    print("Original nonlinear system:")
    system.print()
    poly_sys = polynomialize(system)
    poly_sys.print_substitution_equations()
    poly_sys.print()
    print("=" * 50)
    quadr_system = polynomialize_and_quadratize(system, {u: 1})
    print(quadr_system)