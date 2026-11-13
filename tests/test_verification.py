import pytest
from qbee import *


def heart_model():
    """A pump driven by a phase oscillator (c, s) at rate f = f0 + kN N + g C, with strength u + e C, emptying q."""
    N, C, u, c, s, q = functions("N, C, u, c, s, q")
    d, k = functions("d, k")                       # inputs
    f0, kN, g, e, w, a, b = parameters("f0, kN, g, e, w, a, b")
    f = f0 + kN * N + g * C
    sig = u + e * C
    system = [(N, -N + a * q), (C, d - C), (u, k), (c, -w * f * s), (s, w * f * c), (q, b - sig * f * q * (1 - c))]
    combined = [f * c, f * s, sig * f, sig * f * c, sig * c]
    # quadratize(system, input_free=True), optimal over monomials (about 2 minutes)
    monomial = [C ** 2, c * u, C * c, C * u, N * u, C * N, C ** 2 * c, C * s, C * c * u, N * c * u, N * c, N * s,
                C * N * c]
    return system, combined, monomial


def test_monomial_new_variables():
    x, y = functions("x, y")
    system = [(x, y ** 3), (y, x ** 3)]
    assert is_quadratization(system, [x ** 2, y ** 2, x * y])
    assert not is_quadratization(system, [x ** 2, y ** 2])
    assert reduce_quadratization(system, [x ** 2, y ** 2, x * y, x ** 2 + y ** 2]) == [x ** 2, y ** 2, x * y]


def test_quadratize_result_passes_the_check():
    x, y = functions("x, y")
    system = [(x, y ** 3), (y, x ** 3)]
    res = quadratize(system)
    assert is_quadratization(system, result_new_variables(res, system))


def test_polynomial_new_variables_beat_monomials():
    system, combined, monomial = heart_model()
    assert is_quadratization(system, combined)
    assert is_quadratization(system, monomial)
    assert len(reduce_quadratization(system, combined)) == 5
    for drop in range(5):
        assert not is_quadratization(system, combined[:drop] + combined[drop + 1:])


def test_failures_name_the_equations():
    system, combined, _ = heart_model()
    check = QuadratizationCheck(system, combined)
    q = system[-1][0]
    assert q in check.failures(keep=[0, 1])        # f c, f s alone leave sig * f * q * c cubic


def test_inputs_are_not_allowed_in_new_variables():
    x, u = functions("x, u")
    with pytest.raises(ValueError):
        QuadratizationCheck([(x, x ** 3 * u)], [x * u])


def test_states_given_as_symbols_are_not_parameters():
    import sympy as sp
    p = sp.Symbol("p")
    assert not is_quadratization([(p, p ** 3)], [])
    assert is_quadratization([(p, p ** 3)], [p ** 2])


def test_parameter_only_in_a_new_variable():
    x = functions("x")
    a = parameters("a")
    assert is_quadratization([(x, x ** 3)], [a * x ** 2])


def test_non_polynomial_input_is_rejected():
    import sympy as sp
    x = functions("x")
    with pytest.raises(ValueError):
        QuadratizationCheck([(x, x ** 3)], [1 / x])
    with pytest.raises(ValueError):
        QuadratizationCheck([(x, sp.exp(x))], [])
