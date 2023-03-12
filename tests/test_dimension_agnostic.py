from qbee import *

def test_basic():
    x, Dx = functions("x Dx")
    system = [(x, x**2 * Dx)]
    new_vars = quadratize_dimension_agnostic(system)
    assert len(new_vars) == 2
    x, xt = sp.symbols("x x_tilda")
    assert set(new_vars) == {x**2, x * xt}

def test_scalar_input():
    x, Dx, u = functions("x Dx u")
    system = [(x, x * u * Dx)]
    new_vars = quadratize_dimension_agnostic(system)
    assert len(new_vars) == 1
    x, u = sp.symbols("x u")
    assert new_vars[0] == x * u

def test_scalar_input_free():
    x, Dx, u = functions("x Dx u")
    a = parameters("a")
    system = [(x, x**2 * Dx + a * u * x)]
    new_vars = quadratize_dimension_agnostic(system, input_free=True)
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
