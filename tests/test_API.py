import pytest
import sympy as sp
from qbee import *
from qbee.quadratization import polynomialize_and_quadratize_ode
from qbee.util import derivatives
from functools import partial


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
    res = polynomialize_and_quadratize_ode([
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
    res = polynomialize_and_quadratize_ode(system, {u: 1})
    assert res is not None


def test_polynomialize_and_quadratize_on_already_polynomial_system():
    x, y = functions("x, y")
    res = polynomialize_and_quadratize_ode([
        (x, y ** 5),
        (y, x ** 5)
    ])
    assert res is not None


def test_already_quadratized_system():
    x, y = functions("x, y")
    res = polynomialize_and_quadratize_ode([
        (x, y ** 2),
        (y, x ** 2)
    ])
    assert res.introduced_vars == 0


def test_scalar_pde():
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (u, p * u.diff(x) ** 3 + c),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize(system, pruning_functions=[no_quad_pruning, *default_pruning_rules])
    assert res is not None
    print(res)


def test_scalar_pde_with_elem_funcs():
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (x, x),
        (u, p * u.diff(x) ** 2 * sp.sin(x) + c),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize(system, pruning_functions=[no_quad_pruning, *default_pruning_rules])
    assert res is not None
    print(res)


@pytest.mark.xfail
def test_scalar_functional_pde_without_needed_space_var_derivative():
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (u, p * u.diff(x) ** 2 * sp.sin(x) + c),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    polynomialize_and_quadratize(system, pruning_functions=[no_quad_pruning, *default_pruning_rules])
