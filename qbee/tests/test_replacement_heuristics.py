import pytest
import sympy as sp

from qbee import *
from replacement_heuristics import *
from replacement_heuristics import _compute_monomials_affected, \
    _compute_replacement_value_for_all_monomials, \
    _compute_auxiliary_equation_ql_discrepancy

x, y, z = sp.symbols('x, y, z')
x_dot, y_dot, z_dot = derivatives('x, y, z')


def test_free_variables_count_sorted():
    system = EquationSystem([
        sp.Eq(x_dot, x * y ** 2)
    ])

    replacements = tuple(map(sp.Poly.as_expr, free_variables_count_sorted(system)))
    assert replacements == (y ** 2, x * y)


def test_compute_monomials_affected():
    system = EquationSystem([
        sp.Eq(x_dot, x ** 2 + x ** 2 * y ** 2),
        sp.Eq(y_dot, x * y ** 2 + x ** 2 * y ** 2)
    ])
    replacement = sp.Poly(x * y)

    assert _compute_monomials_affected(system, replacement, unique=False) == 3
    assert _compute_monomials_affected(system, replacement, unique=True) == 2

    replacement = sp.Poly(x ** 3 * y ** 3)
    assert _compute_monomials_affected(system, replacement) == 0


def test_compute_replacement_value_for_all_monomials():
    system = EquationSystem([
        sp.Eq(x_dot, x ** 2 + x ** 2 * y ** 2),
        sp.Eq(y_dot, x * y ** 2 + x ** 2 * y ** 2)
    ])
    replacement = sp.Poly(x * y)

    assert _compute_replacement_value_for_all_monomials(system, replacement) == 3

    replacement = sp.Poly(x ** 3 * y ** 3)
    assert _compute_monomials_affected(system, replacement) == 0


def test_compute_auxiliary_equation_ql_discrepancy():
    system = EquationSystem([
        sp.Eq(x_dot, x ** 2 + x ** 2 * y ** 2),
        sp.Eq(y_dot, -y)
    ])
    replacement = sp.Poly(x * y)

    assert _compute_auxiliary_equation_ql_discrepancy(system, replacement) == 3
