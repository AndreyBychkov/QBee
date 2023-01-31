import pytest
import sympy as sp

from qbee import *
from qbee.polynomialization import eq_list_to_eq_system


def assert_check_poly(expected_system: EquationSystem, actual_system: EquationSystem):
    try:
        assert actual_system.equations == expected_system.equations
    except AssertionError as e:
        if actual_system.is_polynomial():
            raise AssertionError("Systems are not equal but actual system is polynomial.")
        else:
            raise e


def test_already_polynomial():
    x = functions("x")
    p = parameters("p")
    system = eq_list_to_eq_system([
        (x, p * x ** 2 + x ** 3 + 1)
    ])

    res = polynomialize(system)
    assert system.is_polynomial()
    assert system.polynomial_equations == res.polynomial_equations
    assert system.equations == res.equations


def test_sigmoid():
    x = functions("x")
    system = eq_list_to_eq_system([
        (x, 1 / (1 + sp.exp(x)))
    ])

    res = polynomialize(system, keep_laurent=False)
    assert len(res) == 3


def test_sigmoid_inv_arg():
    x = functions("x")
    system = [
        (x, 1 / (1 + sp.exp(1 / x)))
    ]

    res = polynomialize(system, keep_laurent=False)
    assert len(res) == 4


def test_nested_functions():
    x = functions("x")
    system = eq_list_to_eq_system([
        (x, sp.sin(sp.exp(x)))
    ])

    res = polynomialize(system, keep_laurent=False)
    assert len(res) == 4


def test_parameter():
    x, y = functions("x, y")
    p = parameters("p")
    system = [
        (x, sp.exp(p * y)),
        (y, sp.exp(p * x))
    ]
    res = polynomialize(system, upper_bound=4, keep_laurent=False)
    assert all([eq.rhs.is_Mul and sp.Symbol("p") in eq.rhs.args
                for eq in res.polynomial_equations[2:]])


def test_combustion():
    c1, c2, c3, c4, T = functions("c1, c2, c3, c4, T")
    A, Ea, Ru = parameters("A, Ea, Ru")
    eq1 = -A * sp.exp(-Ea / (Ru * T)) * c1 ** 0.2 * c2 ** 1.3
    system = [
        (c1, eq1),
        (c2, 2 * eq1),
        (c3, -eq1),
        (c4, -2 * eq1)
    ]

    res = polynomialize(system, upper_bound=8, keep_laurent=False)
    assert len(res) == 10
