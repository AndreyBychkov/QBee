import pytest
import sympy as sp

from qbee.AST_walk import find_non_polynomial


def test_AST_walk_polynomial_expr():
    x, y = sp.symbols(["x", "y"])

    expr = x ** 2 + y - x * y + y ** 3
    found = find_non_polynomial(expr)
    assert found is None


def test_AST_walk_rational_expr():
    x, y = sp.symbols(["x", "y"])

    expr = (x ** 2 + y) / (x + y ** 2)
    found = find_non_polynomial(expr)
    assert found == 1 / (x + y ** 2)


def test_AST_walk_general_expr():
    x, y = sp.symbols(["x", "y"])

    expr = x ** 2 - y + (x + y) * sp.sin(x)
    found = find_non_polynomial(expr)
    assert found == sp.sin(x)
