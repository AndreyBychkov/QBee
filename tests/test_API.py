import pytest
import sympy as sp
from qbee import *


def test_polynomialize_list_input():
    x, y = sp.symbols("x, y")
    res = polynomialize([
        (x, sp.sin(x) * sp.cos(y)),
        (y, sp.cos(x) * sp.tan(y))
    ])
    assert res is not None


def test_polynomialize_EquationSystem_input():
    x, y = sp.symbols("x, y")
    dx, dy = derivatives([x, y])
    system = EquationSystem([
        sp.Eq(dx, sp.sin(x) * sp.cos(y)),
        sp.Eq(dy, sp.cos(x) * sp.tan(y))
    ])
    res = polynomialize(system)
    assert res is not None


def test_polynomialize_and_quadratize_list_input():
    x, y = sp.symbols("x, y")
    res = polynomialize_and_quadratize([
        (x, sp.sin(x) * sp.cos(y)),
        (y, sp.cos(x) * sp.tan(y))
    ])
    assert res is not None


def test_polynomialize_and_quadratize_EquationSystem_input():
    x, y = sp.symbols("x, y")
    dx, dy = derivatives([x, y])
    system = EquationSystem([
        sp.Eq(dx, sp.sin(x) * sp.cos(y)),
        sp.Eq(dy, sp.cos(x) * sp.tan(y))
    ])
    res = polynomialize_and_quadratize(system)
    assert res is not None


def test_polynomialize_and_quadratize_on_already_polynomial_system():
    x, y = sp.symbols("x, y")
    res = polynomialize_and_quadratize([
        (x, y**5),
        (y, x**5)
    ])
    assert res is not None
