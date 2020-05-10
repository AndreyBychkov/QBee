import pytest

from pytest import raises
from structures import *

x = sp.symbols('x')
x_dot = make_derivative_symbol(x)


def test_monomial_extraction():
    system = EquationSystem([sp.Eq(x_dot, 1 + x ** 2 + x ** 3)])
    assert system.monomials == (1, x ** 2, x ** 3)

    system = EquationSystem([
        sp.Eq(x_dot, 1),
        sp.Eq(x_dot, x),
        sp.Eq(x_dot, x**2 - x**3)
    ])
    assert system.monomials == (1, x, x**2, -x**3)

    with raises(AssertionError):
        system = EquationSystem([sp.Eq(x_dot, 1 / x)])
        system.monomials
