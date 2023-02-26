from qbee import *
from qbee import polynomialize_and_quadratize
from qbee.polynomialization import eq_list_to_eq_system


def test_handmade_negative_laurent():
    x = functions("x")
    system = eq_list_to_eq_system([
        (x, 1 / x ** 2)
    ])

    assert system.is_polynomial()


def test_sigmoid_inv_arg():
    x = functions("x")
    system = [
        (x, 1 / (1 + sp.exp(1 / x)))
    ]

    res = polynomialize_and_quadratize(system)
    assert res.new_vars_count == 5


def test_solarwind_quad_only():
    coef_field = sp.FractionField(sp.QQ, ["c"])
    R, v1, v2, v3 = sp.ring("v1, v2, v3", coef_field)
    c, *_ = coef_field.gens
    # vt = vx * vinv - c * vx
    v1t = (v1 + v2) * v1 ** (-1) - c * (v1 + v2)
    v2t = (v1 + v2 + v3) * v2 ** (-1) + c * (v1 + v2 + v3)
    v3t = (v2 + v3) * v3 ** (-1) + c * (v2 + v3)

    quadr_system = quadratize([v1t, v2t, v3t])
    if quadr_system:
        print(quadr_system)
    assert quadr_system.new_vars_count == 11


def test_solarwind():
    v1, v2, v3 = functions("v1, v2, v3")
    c = parameters("c")
    # vt = vx * vinv - c * vx
    v1t = (v1 + v2) * v1 ** (-1) - c * (v1 + v2)
    v2t = (v1 + v2 + v3) * v2 ** (-1) + c * (v1 + v2 + v3)
    v3t = (v2 + v3) * v3 ** (-1) + c * (v2 + v3)
    system = [
        (v1, v1t),
        (v2, v2t),
        (v3, v3t)
    ]

    quadr_system = polynomialize_and_quadratize(system)
    if quadr_system:
        print(quadr_system)
    assert quadr_system.new_vars_count == 11
