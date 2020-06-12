import pytest
import sympy as sp

from qbee import polynomialize, quadratic_linearize, EquationSystem, derivatives

x, y, z = sp.symbols('x, y, z')
dot_x, dot_y, dot_z = derivatives('x, y, z')


def test_sigmoid():
    system = EquationSystem([
        sp.Eq(dot_x, 1 / (1 + sp.exp(x)))
    ])
    poly_system = polynomialize(system)
    ql_result = quadratic_linearize(poly_system)
    assert len(ql_result.system.equations) == 4


def test_zero_system():
    w = sp.symbols('w')
    dot_w = derivatives('w')

    system = EquationSystem([
        sp.Eq(dot_x, 0),
        sp.Eq(dot_y, 0),
        sp.Eq(dot_z, 0),
        sp.Eq(dot_w, x ** 2 * y ** 2 * z ** 2)
    ])

    poly_system = polynomialize(system)
    ql_result = quadratic_linearize(poly_system)

    assert len(ql_result.system.equations) == 5


def test_x_sigmoid():
    system = EquationSystem([
        sp.Eq(dot_x, x / (1 + sp.exp(x)))
    ])
    poly_system = polynomialize(system)
    ql_result = quadratic_linearize(poly_system, limit_depth=2, initial_max_depth=2)
    assert len(ql_result.system.equations) == 5


def test_rabinovich_fabrikant():
    a, b = sp.symbols('a, b')

    system = EquationSystem([
        sp.Eq(dot_x, y * (z - 1 + x ** 2) + a * x),
        sp.Eq(dot_y, x * (3 * z + 1 - x ** 2) + a * y),
        sp.Eq(dot_z, -2 * z * (b + x * y))
    ], parameter_variables=[a, b])

    ql_res = quadratic_linearize(system, initial_max_depth=3, limit_depth=3, heuristics='auxiliary-equation-ql-discrepancy')
    assert len(ql_res.system.equations) == 6
