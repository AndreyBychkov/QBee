"""
Checking and reducing a quadratization with polynomial new variables (qbee.verification).

The system is a pump driven by a phase oscillator (c, s) at rate f = f0 + kN N + g C, with strength u + e C,
emptying a volume q; d and k are inputs. The optimal quadratization over monomials has 13 new variables
(quadratize, about 2 minutes). Products of the model's own expressions need 5:
f c, f s, (u + e C) f, (u + e C) f c, (u + e C) c.
quadratize_affine finds such a set automatically: it makes the affine sub-expressions f and u + e C
coordinates and runs the same search there (about 1 second).
"""
import time

from qbee import *

if __name__ == '__main__':
    N, C, u, c, s, q = functions("N, C, u, c, s, q")
    d, k = functions("d, k")
    f0, kN, g, e, w, a, b = parameters("f0, kN, g, e, w, a, b")
    f = f0 + kN * N + g * C
    sig = u + e * C
    system = [(N, -N + a * q), (C, d - C), (u, k), (c, -w * f * s), (s, w * f * c), (q, b - sig * f * q * (1 - c))]

    combined = [f * c, f * s, sig * f, sig * f * c, sig * c]
    print("polynomial new variables:", len(combined), "quadratize:", is_quadratization(system, combined))
    print("irredundant:", len(reduce_quadratization(system, combined)))

    t0 = time.time()
    res = quadratize(system, input_free=True)
    monomial = result_new_variables(res, system)
    print(f"optimal monomial quadratization: {len(monomial)} new variables ({time.time() - t0:.0f} s),",
          "checked:", is_quadratization(system, monomial))

    t0 = time.time()
    auto = quadratize_affine(system)
    print(f"automatic (affine sub-expressions as coordinates): {len(auto.new_vars)} new variables "
          f"({time.time() - t0:.0f} s), checked:", is_quadratization(system, auto.new_vars))
    for y, form in auto.coordinates.items():
        print(f"    coordinate {y} = {form}")
