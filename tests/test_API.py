import pytest

import sympy as sp
from qbee import *
from qbee.polynomialization import EquationSystem
from qbee.experimental import polynomialize_and_quadratize_pde
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
    system = EquationSystem({
        x: sp.sin(k * x + u),
        y: p * sp.cos(y)
    }, [p, k], [u])
    res = polynomialize(system)
    assert len(res) > 2


def test_polynomialize_and_quadratize_list_input():
    x, y, u = functions("x, y, u", real=True)  # Identical to sympy.symbols
    p, k = parameters("p, k")
    res = polynomialize_and_quadratize([
        (x, sp.sin(k * x + u) * y),
        (y, p * sp.cos(y))
    ], input_der_orders={u: 1})
    assert res is not None


def test_polynomialize_and_quadratize_EquationSystem_input():
    x, y, u = sp.symbols("x, y, u")
    p, k = sp.symbols("p, k")
    system = EquationSystem({
        x: sp.sin(k * x + u) * y,
        y: p * sp.cos(y)
    }, [p, k], [u])
    res = polynomialize_and_quadratize(system, input_der_orders={u: 1})
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
    assert res.new_vars_count == 0


def test_pruning_by_quadratizations_found():
    x, y = functions("x, y")
    system = [(x, y ** 3), (y, x ** 3)]
    optimal = quadratize(system)
    first = quadratize(system, pruning_functions=[pruning_by_quadratizations_found, *default_pruning_rules])
    assert first.quadratization.is_quadratized()
    assert first.nodes_traversed < optimal.nodes_traversed
    assert len(first.introduced_variables) >= len(optimal.introduced_variables)

    many = quadratize(system, pruning_functions=[partial(pruning_by_quadratizations_found, quadratizations=10 ** 6),
                                                 *default_pruning_rules])
    assert len(many.introduced_variables) == len(optimal.introduced_variables)


def test_first_quadratization_with_inputs():
    from examples import generate_cardiovascular_respiratory
    system = generate_cardiovascular_respiratory(gas_exchange=False)   # 11 states, 5 inputs; about 4 s
    res = quadratize(system, input_free=True, calc_upper_bound=False,
                     pruning_functions=[pruning_by_quadratizations_found, *default_pruning_rules])
    assert res is not None and res.quadratization.is_quadratized()
    assert res.nodes_traversed < 200                                    # 87: the search stopped early


def test_first_quadratization_with_upper_bound():
    # the upper-bound search's quadratization bounds the main search; one that does not beat the bound must not count
    from examples import generate_cubic_cycle
    res = quadratize(generate_cubic_cycle(4), calc_upper_bound=True,
                     pruning_functions=[pruning_by_quadratizations_found, *default_pruning_rules])
    assert res is not None and res.quadratization.is_quadratized()


def test_quadratizations_found_resets_between_searches():
    from qbee.quadratization import BranchAndBound, PolynomialSystem, pruning_by_best_nvars
    R, a, b = sp.ring("a, b", sp.QQ)
    algo = BranchAndBound(PolynomialSystem([b ** 3, a ** 3]), None, default_generation, default_scoring,
                          [pruning_by_best_nvars, pruning_by_quadratizations_found])
    first = algo.quadratize().num_introduced_vars
    assert algo.quadratize().num_introduced_vars == first
