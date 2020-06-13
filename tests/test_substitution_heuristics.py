import pytest
import sympy as sp

from qbee import EquationSystem, derivatives
from qbee.substitution_heuristics import free_variables_count_sorted, \
    _compute_monomials_affected, \
    _compute_substitution_value_for_all_monomials, \
    _compute_auxiliary_equation_quadratic_discrepancy

x, y, z = sp.symbols('x, y, z')
x_dot, y_dot, z_dot = derivatives('x, y, z')


def test_free_variables_count_sorted():
    system = EquationSystem([
        sp.Eq(x_dot, x * y ** 2)
    ])

    substitutions = tuple(map(sp.Poly.as_expr, free_variables_count_sorted(system)))
    assert substitutions == (y ** 2, x * y)


def test_compute_monomials_affected():
    system = EquationSystem([
        sp.Eq(x_dot, x ** 2 + x ** 2 * y ** 2),
        sp.Eq(y_dot, x * y ** 2 + x ** 2 * y ** 2)
    ])
    substitution = sp.Poly(x * y)

    assert _compute_monomials_affected(system, substitution, unique=False) == 3
    assert _compute_monomials_affected(system, substitution, unique=True) == 2

    substitution = sp.Poly(x ** 3 * y ** 3)
    assert _compute_monomials_affected(system, substitution) == 0


def test_compute_replacement_value_for_all_monomials():
    system = EquationSystem([
        sp.Eq(x_dot, x ** 2 + x ** 2 * y ** 2),
        sp.Eq(y_dot, x * y ** 2 + x ** 2 * y ** 2)
    ])
    substitution = sp.Poly(x * y)

    assert _compute_substitution_value_for_all_monomials(system, substitution) == 3

    substitution = sp.Poly(x ** 3 * y ** 3)
    assert _compute_monomials_affected(system, substitution) == 0


def test_compute_auxiliary_equation_ql_discrepancy():
    system = EquationSystem([
        sp.Eq(x_dot, x ** 2 + x ** 2 * y ** 2),
        sp.Eq(y_dot, -y)
    ])
    substitution = sp.Poly(x * y)

    assert _compute_auxiliary_equation_quadratic_discrepancy(system, substitution) == 3
