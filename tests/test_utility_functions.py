import pytest
import sympy as sp

from qbee import derivatives
from qbee.util import is_monomial_divisor, symbol_from_derivative

x, y, z = sp.symbols('x, y, z')


def test_is_monomial_divisor():
    numerator = x ** 2 * y ** 3
    denominator = x * y
    assert is_monomial_divisor(numerator, denominator) is True

    numerator = x ** 2 * y ** 3
    denominator = x * y ** 4
    assert is_monomial_divisor(numerator, denominator) is False


def test_symbol_from_derivative():
    y0 = sp.Symbol("y_{0}")
    dy0 = derivatives(y0)

    actual = symbol_from_derivative(dy0)
    expected = y0

    assert actual == expected
