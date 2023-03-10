from qbee import *

def test_basic():
    x, Dx = functions("x Dx")
    system = [(x, x**2 * Dx)]
    new_vars = quadratize_dimension_agnostic(system)
    assert len(new_vars) == 2
    x, xt = sp.symbols("x x_tilda")
    assert set(new_vars) == {x**2, x * xt}

def test_solarwind():
    u, v, Du, Dv = functions("u v Du Dv")
    C1 = parameters("C1")
    system = [
        (v, Du + C1 * Dv),
        (u, Du / v + C1 * Dv / v)
    ]
    new_vars = quadratize_dimension_agnostic(system)
    assert len(new_vars) == 4
    u, v, ut, vt = sp.symbols("u v u_tilda v_tilda")
    assert set(new_vars) == {1 / v, vt / v, u / v, ut / v}
