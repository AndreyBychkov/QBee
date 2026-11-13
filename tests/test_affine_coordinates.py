import sympy as sp

from qbee import functions, parameters
from qbee.affine_coordinates import affine_subexpressions, choose_coordinates, quadratize_affine, transform
from qbee.verification import is_quadratization


def pump():
    """A pump driven by a phase oscillator (c, s) at rate f, with strength u + e C, emptying a volume q."""
    N, C, u, c, s, q = functions("N, C, u, c, s, q")
    d, k = functions("d, k")
    f0, kN, g, e, w, a, b = parameters("f0, kN, g, e, w, a, b")
    f = f0 + kN * N + g * C
    sig = u + e * C
    system = [(N, -N + a * q), (C, d - C), (u, k), (c, -w * f * s), (s, w * f * c), (q, b - sig * f * q * (1 - c))]
    return system, f, sig


def test_affine_subexpressions_by_frequency_without_inputs():
    system, f, sig = pump()
    found = dict(affine_subexpressions(system))
    assert found[f] == 3                                    # in c', s' and q'
    assert sig in found
    C, d = functions("C, d")
    assert not any(e.has(d) for e in found)                 # C' = d - C contains the input d
    assert all(len(e.free_symbols & {x for x, _ in system}) + len(e.atoms(sp.Function)) > 0 for e in found)


def test_dependent_forms_are_not_both_coordinates():
    N, C = functions("N, C")
    system = [(N, -(N + C) * (2 * N + 2 * C + 1)), (C, -C)]
    chosen = choose_coordinates(system, candidates=[N + C, 2 * N + 2 * C + 1])
    assert len(chosen) == 1


def test_transform_is_the_chain_rule():
    system, f, sig = pump()
    chosen = choose_coordinates(system)
    new_system, coords, replaced = transform(system, chosen)
    back = dict(coords)
    rates = {x: sp.sympify(r) for x, r in system}
    states = [x for x, _ in system]
    for (x, _), (y, r) in zip(system, new_system):
        if y in coords:
            want = sum(sp.diff(coords[y], z) * rates[z] for z in states)
        else:
            want = rates[x]
        assert sp.expand(r.xreplace(back) - want) == 0, x


def test_pump_needs_five_polynomial_new_variables():
    """QBee's optimal monomial quadratization of the pump has 13 new variables; in coordinates where f and u + e C
    replace states it has 5, which mapped back quadratize the original system."""
    system, f, sig = pump()
    res = quadratize_affine(system)
    assert len(res.new_vars) == 5
    assert is_quadratization(system, res.new_vars)
    assert f in res.coordinates.values() and sig in res.coordinates.values()
