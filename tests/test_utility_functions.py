import pytest
import sympy as sp

from util import *

x, y, z = sp.symbols('x, y, z')


def test_is_monomial_divisor():
    numerator = x ** 2 * y ** 3
    denominator = x * y
    assert is_monomial_divisor(numerator, denominator) is True

    numerator = x ** 2 * y ** 3
    denominator = x * y ** 4
    assert is_monomial_divisor(numerator, denominator) is False
