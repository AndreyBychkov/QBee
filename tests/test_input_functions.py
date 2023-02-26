import numpy as np
import pytest
from examples import *
from functools import partial
from qbee.experimental import to_odeint


@pytest.mark.xfail
def test_1d():
    x, u = functions("x, u")
    system = [
        (x, x ** 2 * u),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=10000)
    res = quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)


@pytest.mark.xfail
def test_2d_xy_in_x():
    x, y, u = functions("x, y, u")
    system = [
        (x, x * y * u),
        (y, y ** 2),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=20000)
    res = quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)


def test_2d_y2_in_x_and_x2_in_y():
    x, y, u = functions("x, y, u")
    system = [
        (x, y ** 2 * u),
        (y, x ** 2),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=20000)
    res = quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)

@pytest.mark.expensive
def test_2d_y2_in_x_and_x4_in_y():
    x, y, u = functions("x, y, u")
    system = [
        (x, y ** 2 * u),
        (y, x ** 4),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=20000)
    res = quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)


@pytest.mark.expensive
def test_2d_y2_in_x_and_x4_in_y_with_z_eq_y2():
    x, y, z, u = functions("x, y, z, u")
    system = [
        (x, z * u),
        (y, x ** 4),
        (z, 2 * x ** 4 * y ** 2)
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=20000)
    res = quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)


@pytest.mark.experimental
def test_numerical():
    x, u = functions("x, u")
    p = parameters("p")
    system = [
        (x, p * sp.sin(x) + u),
    ]
    res = polynomialize_and_quadratize(system, {u: 0})
    print(res)

    ts = np.linspace(0, 100, 10000)
    t = INDEPENDENT_VARIABLE
    my_odeint = to_odeint(system, {"x": 0.1}, {"p": 0.1}, {"u": sp.sin(t)})
    my_odeint_quad = to_odeint(res, {"x": 0.1}, {"p": 0.1}, {"u": sp.sin(t)})

    num_res = my_odeint(ts, atol=1e-11)
    num_res_quad = my_odeint_quad(ts, atol=1e-11)

    assert np.linalg.norm(num_res[:, 0] - num_res_quad[:, 0], np.inf) < 1e-6
