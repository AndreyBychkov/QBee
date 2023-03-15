import pytest

from functools import partial
from qbee import *
from qbee.experimental import polynomialize_and_quadratize_pde, multivariable_functions


@pytest.mark.experimental
def test_scalar_pde():
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (x, x),
        (u, p * u.diff(x) ** 3 + c),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize_pde(system, pruning_functions=[no_quad_pruning, *default_pruning_rules])
    assert res is not None
    print(res)


@pytest.mark.experimental
def test_scalar_pde_with_elem_funcs():
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (x, x),
        (u, p * u.diff(x) ** 2 * sp.sin(x) + c),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize_pde(system, pruning_functions=[no_quad_pruning, *default_pruning_rules])
    assert res is not None
    print(res)


@pytest.mark.experimental
@pytest.mark.xfail
def test_scalar_functional_pde_without_needed_space_var_derivative():
    x, c = functions("x, c")
    u = multivariable_functions("u", [x])
    p = parameters("p")
    system = [
        (u, p * u.diff(x) ** 2 * sp.sin(x) + c),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    polynomialize_and_quadratize_pde(system, pruning_functions=[no_quad_pruning, *default_pruning_rules])


@pytest.mark.experimental
def test_dym_equations():
    x = functions("x")
    u = multivariable_functions("u", [x])
    rhs = u ** 3 * u.diff(x, 3)
    system = [
        (x, x),
        (u, rhs),
        (u.diff(x), rhs.diff(x)),
        (u.diff(x, 2), rhs.diff(x, 2))
    ]
    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=1000)
    res = polynomialize_and_quadratize_pde(system, pruning_functions=[*default_pruning_rules,
                                                                      no_quad_pruning
                                                                      ])
    assert res is not None
    print(res)
