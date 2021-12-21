import itertools
import sympy as sp
from qbee import *
from functools import reduce
from random import randrange
from operator import add, mul


def generate_x_plus_1(deg):
    x = functions("x")
    return [(x, (x + 1) ** deg)]


def generate_dense_polynomial(nvars, deg):
    R, *variables = sp.ring([f"x{i}" for i in range(nvars)], sp.QQ)

    def gen_eq_deg(d: int):
        return reduce(add, [randrange(-10, 10) * reduce(mul, p) for p in itertools.product(variables, repeat=d)]) \
               + randrange(-10, 10)

    def gen_eq():
        return reduce(add, [gen_eq_deg(d) for d in range(1, deg + 1)])

    return [gen_eq() for _ in range(nvars)]


def generate_circular(deg):
    x, y = functions("x, y")
    return [(x, y ** deg), (y, x ** deg)]


def generate_hard(deg):
    a, b, c = functions("a, b, c")
    return [(a, c ** deg + a ** 2 * b ** 2 * c ** 3),
            (b, a ** 2),
            (c, b ** 2)]


def generate_hill(k):
    h, i, t = functions("h, i, t")
    return [(h, k * i ** 2 * t ** (k - 1)),
            (i, -k * i ** 2 * t ** (k - 1)),
            (t, sp.sympify(1))]


def generate_lifeware_conjecture(n):
    """
    Generate the system from Conjecture 3.1 from https://hal.inria.fr/hal-02900798v2/document
    Also is called Long Monomial in benchmarks.
    """
    variables = functions(", ".join(["x" + str(i) for i in range(1, n + 1)]))
    prod_all = 1
    for i in range(n):
        prod_all *= variables[i] ** 2

    system = []
    for i in range(n):
        system.append((variables[i], variables[(i + 1) % n] ** 2 + prod_all))
    return system


def generate_cubic_cycle(n):
    variables = functions(", ".join(["x" + str(i) for i in range(1, n + 1)]))
    system = []
    for i in range(n):
        system.append((variables[i], variables[(i + 1) % n] ** 3))
    return system


def generate_cubic_bicycle(n):
    variables = functions(", ".join(["x" + str(i) for i in range(1, n + 1)]))
    system = []
    for i in range(n):
        system.append((variables[i], variables[(i + 1) % n] ** 3 + variables[(i - 1) % n] ** 3))
    return system


def generate_selkov(a, b):
    x, y = functions("x, y")
    return [(x, -x + a * y + x ** 2 * y),
            (y, b - a * y - x ** 2 * y)]


def generate_exponential():
    coef_field = sp.FractionField(sp.QQ, ["Apsi", "Atheta", "bpsi", "btheta", "gamma", "B", "D"])
    Apsi, Atheta, bpsi, btheta, gamma, B, D = coef_field.gens
    R, psi, theta, exp, thinv2, u, du = sp.ring(["psi", "theta", "exp", "thinv2", "u", "u'"], coef_field)
    eqs = [
        Apsi * psi + bpsi * u - D * psi * exp,
        Atheta * theta + btheta * u + B * D * psi * exp
    ]
    eqs.append(gamma * eqs[0] * thinv2 * exp)
    eqs.append(eqs[0] * (-2) * thinv2 * psi)
    eqs.append(du)
    eqs.append(du ** 2)
    return eqs

for i in range(6, 8):
    sys = generate_cubic_cycle(i)
    quadr_system = polynomialize_and_quadratize(sys, offset=1)
    print(quadr_system)
