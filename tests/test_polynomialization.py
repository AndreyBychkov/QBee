import pytest
import sympy as sp

from qbee import EquationSystem, derivatives, polynomialize

x, y, z = sp.symbols('x, y, z')
dot_x, dot_y, dot_z = derivatives('x, y, z')


def assert_check_poly(expected_system: EquationSystem, actual_system: EquationSystem):
    try:
        assert actual_system.equations == expected_system.equations
    except AssertionError as e:
        if actual_system.is_polynomial("full"):
            raise AssertionError("Systems are not equal but actual system is polynomial.")
        else:
            raise e


def test_already_polynomial():
    system = EquationSystem([
        sp.Eq(dot_x, x + x ** 2 + 3)
    ])

    assert polynomialize(system).equations == system.equations


def test_sigmoid_diff():
    system = EquationSystem([
        sp.Eq(dot_x, 1 / (1 + sp.exp(x)))
    ])

    poly_system = polynomialize(system)

    _, y0, y1 = poly_system.variables.variables
    dot_y0, dot_y1 = derivatives([y0, y1])
    expected_system = EquationSystem([
        sp.Eq(dot_x, y1),
        sp.Eq(dot_y0, y0 * y1),
        sp.Eq(dot_y1, -y0 * y1 ** 3)
    ])

    assert_check_poly(expected_system, poly_system)


def test_parameter():
    k = sp.Symbol('k')

    system = EquationSystem([
        sp.Eq(dot_x, sp.exp(k * x)),
        sp.Eq(dot_y, sp.exp(k * x))
    ], parameter_variables=[k])

    poly_system = polynomialize(system)

    _, y0 = poly_system.variables.variables
    dot_y0 = derivatives(y0)
    expected_system = EquationSystem([
        sp.Eq(dot_x, y0),
        sp.Eq(dot_y, y0),
        sp.Eq(dot_y0, k * y0 ** 2)
    ], parameter_variables=[k])

    assert_check_poly(expected_system, poly_system)
