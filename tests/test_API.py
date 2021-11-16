import pytest
import sympy as sp
from qbee import *
from qbee.util import derivatives


def test_polynomialize_list_input():
    x, y, u = functions("x, y, u", real=True)
    p, k = parameters("p, k")
    res = polynomialize([
        (x, sp.sin(k * x + u)),
        (y, p * sp.cos(y))
    ])
    assert len(res) > 2


def test_polynomialize_EquationSystem_input():
    x, y, u = sp.symbols("x, y, u")
    p, k = sp.symbols("p, k")
    dx, dy = derivatives([x, y])
    system = EquationSystem([
        sp.Eq(dx, sp.sin(k * x + u)),
        sp.Eq(dy, p * sp.cos(y))
    ], [p, k], [u])
    res = polynomialize(system)
    assert len(res) > 2


def test_polynomialize_and_quadratize_list_input():
    x, y, u = functions("x, y, u", real=True)  # Identical to sympy.symbols
    p, k = parameters("p, k")
    res = polynomialize_and_quadratize([
        (x, sp.sin(k * x + u) * y),
        (y, p * sp.cos(y))
    ], {u: 1})
    assert res is not None


def test_polynomialize_and_quadratize_EquationSystem_input():
    x, y, u = sp.symbols("x, y, u")
    p, k = sp.symbols("p, k")
    dx, dy = derivatives([x, y])
    system = EquationSystem([
        sp.Eq(dx, sp.sin(k * x + u) * y),
        sp.Eq(dy, p * sp.cos(y))
    ], [p, k], [u])
    res = polynomialize_and_quadratize(system, {u: 1})
    assert res is not None


def test_polynomialize_and_quadratize_on_already_polynomial_system():
    x, y = functions("x, y")
    res = polynomialize_and_quadratize([
        (x, y ** 5),
        (y, x ** 5)
    ])
    assert res is not None


def test_already_quadratized_system():
    x, y = functions("x, y")
    res = polynomialize_and_quadratize([
        (x, y ** 2),
        (y, x ** 2)
    ])
    assert res.introduced_vars == 0

