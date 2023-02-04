import pytest
import sympy as sp
from qbee import *
from qbee.polynomialization import eq_list_to_eq_system


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

    res = polynomialize(system)
    assert len(res) == 3


def test_sigmoid_inv_arg():
    x = functions("x", laurent=False)
    system = [
        (x, 1 / (1 + sp.exp(1 / x)))
    ]

    res = polynomialize(system)
    assert len(res) == 4


def test_nested_functions():
    x = functions("x", laurent=False)
    system = eq_list_to_eq_system([
        (x, sp.sin(sp.exp(x)))
    ])

    res = polynomialize(system)
    assert len(res) == 4


def test_parameter():
    x, y = functions("x, y", laurent=False)
    p = parameters("p")
    system = [
        (x, sp.exp(p * y)),
        (y, sp.exp(p * x))
    ]
    res = polynomialize(system, upper_bound=4)
    assert all([eq.rhs.is_Mul and sp.Symbol("p") in eq.rhs.args
                for eq in res.polynomial_equations[2:]])


def test_combustion():
    c1, c2, c3, c4, T = functions("c1, c2, c3, c4, T", laurent=False)
    A, Ea, Ru = parameters("A, Ea, Ru")
    eq1 = -A * sp.exp(-Ea / (Ru * T)) * c1 ** 0.2 * c2 ** 1.3
    system = [
        (c1, eq1),
        (c2, 2 * eq1),
        (c3, -eq1),
        (c4, -2 * eq1)
    ]

    res = polynomialize(system, upper_bound=8)
    assert len(res) == 10
