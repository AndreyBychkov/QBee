import pytest

from sympy import ring, QQ
from qbee import *


def test_simple():
    R, x, y = ring(['x', 'y'], QQ)
    poly_system = PolynomialSystem([x ** 2 + y, 3 * x * y ** 3 - y])
    algo = BranchAndBound(poly_system, 10, heuristics=default_score)
    res = algo.quadratize()
    assert res.system.is_quadratized()
    assert res.introduced_vars <= 2
    print(res)


def test_poly():
    R, x, y = ring(['x', 'y'], QQ)
    poly_system = PolynomialSystem([(x + 1) ** 8, y])
    algo = BranchAndBound(poly_system, 20, heuristics=default_score)
    res = algo.quadratize()
    assert res.system.is_quadratized()
    assert res.introduced_vars <= 4
    print(res)


def test_poly_long():
    R, x, y = ring(['x', 'y'], QQ)
    poly_system = PolynomialSystem([(x + 1) ** 15, y])
    algo = BranchAndBound(poly_system, 8, heuristics=default_score)
    res = algo.quadratize()
    assert res.system.is_quadratized()
    assert res.introduced_vars <= 7
    print(res)


def test_xSigmoid():
    R, x, y, z = ring(["x", "y", "z"], QQ)
    poly_system = PolynomialSystem([x * z, x * y * z, x * y * z ** 3])
    algo = BranchAndBound(poly_system, 5, heuristics=default_score)
    res = algo.quadratize()
    assert res.system.is_quadratized()
    assert res.introduced_vars <= 2
    print(res)


def test_RabinovichFabricant():
    R, x, y, z = ring(["x", "y", "z"], QQ)
    poly_system = PolynomialSystem([y * (z - 1 - x ** 2) + x, x * (3 * z + 1 - x ** 2) + y, -2 * z * (2 + x * y)])
    algo = BranchAndBound(poly_system, 20, heuristics=default_score)
    res = algo.quadratize()
    assert res.system.is_quadratized()
    assert res.introduced_vars <= 3
    print(res)
