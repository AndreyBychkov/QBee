import pytest
from examples import *
from functools import partial


@pytest.mark.xfail
def test_1d():
    x, u = functions("x, u")
    system = [
        (x, x ** 2 * u),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=10000)
    res = polynomialize_and_quadratize(system, {u: 0}, pruning_functions=[
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
    res = polynomialize_and_quadratize(system, {u: 0}, pruning_functions=[
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
    res = polynomialize_and_quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)


def test_2d_y2_in_x_and_x4_in_y():
    x, y, u = functions("x, y, u")
    system = [
        (x, y ** 2 * u),
        (y, x ** 4),
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=20000)
    res = polynomialize_and_quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)


def test_2d_y2_in_x_and_x4_in_y_with_z_eq_y2():
    x, y, z, u = functions("x, y, z, u")
    system = [
        (x, z * u),
        (y, x ** 4),
        (z, 2 * x ** 4 * y ** 2)
    ]

    no_quad_pruning = partial(pruning_by_nodes_without_quadratization_found, nodes_processed=20000)
    res = polynomialize_and_quadratize(system, {u: 0}, pruning_functions=[
        no_quad_pruning,
        *default_pruning_rules])
    assert res is not None
    print(res)
