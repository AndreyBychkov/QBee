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


def test_handmade_sin_cos():
    x = functions("x")
    system = eq_list_to_eq_system([
        (x, sp.sin(x))
    ])

    xs = sp.Symbol("x")
    system.add_new_var(sp.sin(xs))
    system.add_new_var(sp.cos(xs))
    assert system.is_polynomial()


def test_handmade_sin_cos_inverted_order():
    x = functions("x")
    system = eq_list_to_eq_system([
        (x, sp.sin(x))
    ])

    xs = sp.Symbol("x")
    system.add_new_var(sp.cos(xs))
    system.add_new_var(sp.sin(xs))
    assert system.is_polynomial()  # Should it be correct?


def test_handmade_sigmoid():
    x = functions("x")
    system = eq_list_to_eq_system([
        (x, 1 / (1 + sp.exp(x)))
    ])

    xs = sp.Symbol("x")
    system.add_new_var(sp.exp(xs))
    system.add_new_var(1 / (1 + sp.exp(xs)))
    assert system.is_polynomial()


def test_handmade_negative():
    x = functions("x", laurent=False)
    system = eq_list_to_eq_system([
        (x, 1 / x ** 2)
    ])
    assert not system.is_polynomial()

    xs = sp.Symbol("x")
    system.add_new_var(1 / xs)
    assert system.is_polynomial()


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

def test_parametric_pow():
    x = functions("x")
    a = parameters("a")
    system = [(x, x**a)]
    res = polynomialize(system)
    assert len(res) == 2

    x = functions("x", laurent=False)
    res = polynomialize([(x, x**a)])
    assert len(res) == 3

    res = polynomialize([(x, x**(-a))])
    assert len(res) == 3

def test_hill():
    x = functions("x")
    a, b = parameters("a b")
    system = [(x, 1 / (a + x**b))]
    res = polynomialize(system)
    assert len(res) == 3

    x = functions("x", laurent=False)
    res = polynomialize([(x, 1 / (a + x**b))])
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
