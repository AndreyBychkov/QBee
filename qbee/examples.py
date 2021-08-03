import itertools
from functools import partial, reduce
from random import randrange
from operator import add, mul
from .quadratization import *


def generate_x_plus_1(deg):
    R, x = ring(['x'], QQ)
    return [(x + 1) ** deg]


def generate_dense_polynomial(nvars, deg):
    R, *variables = ring([f"x{i}" for i in range(nvars)], QQ)

    def gen_eq_deg(d: int):
        return reduce(add, [randrange(-10, 10) * reduce(mul, p) for p in itertools.product(variables, repeat=d)]) \
               + randrange(-10, 10)

    def gen_eq():
        return reduce(add, [gen_eq_deg(d) for d in range(1, deg + 1)])

    return [gen_eq() for _ in range(nvars)]


def generate_circular(deg):
    R, x, y = ring(['x', 'y'], QQ)
    return [y ** deg, x ** deg]


def generate_hard(deg):
    R, a, b, c = ring(['a', 'b', 'c'], QQ)
    return [c ** deg + a ** 2 * b ** 2 * c ** 3,
            a ** 2,
            b ** 2]


def generate_hill(k):
    R, h, i, t = ring(['h', 'i', 't'], QQ)
    return [k * i ** 2 * t ** (k - 1),
            -k * i ** 2 * t ** (k - 1),
            R.one]


def generate_lifeware_conjecture(n):
    """
    Generate the system from Conjecture 3.1 from https://hal.inria.fr/hal-02900798v2/document
    Also is called Long Monomial in benchmarks.
    """
    variables = ring([f"x{i}" for i in range(n)], QQ)[1:]
    prod_all = 1
    for i in range(n):
        prod_all *= variables[i] ** 2

    system = []
    for i in range(n):
        system.append(variables[(i + 1) % n] ** 2 + prod_all)
    return system


def generate_cubic_cycle(n):
    variables = ring([f"x{i}" for i in range(n)], QQ)[1:]
    system = []
    for i in range(n):
        system.append(variables[(i + 1) % n] ** 3)
    return system


def generate_cubic_bicycle(n):
    variables = ring([f"x{i}" for i in range(n)], QQ)[1:]
    system = []
    for i in range(n):
        system.append(variables[(i + 1) % n] ** 3 + variables[(i - 1) % n] ** 3)
    return system


def generate_selkov(a, b):
    R, x, y = ring(['x', 'y'], QQ)
    return [-x + a * y + x ** 2 * y,
            b - a * y - x ** 2 * y]


def generate_exponential():
    coef_field = FractionField(QQ, ["Apsi", "Atheta", "bpsi", "btheta", "gamma", "B", "D"])
    Apsi, Atheta, bpsi, btheta, gamma, B, D = coef_field.gens
    R, psi, theta, exp, thinv2, u, du = ring(["psi", "theta", "exp", "thinv2", "u", "u'"], coef_field)
    eqs = [
        Apsi * psi + bpsi * u - D * psi * exp,
        Atheta * theta + btheta * u + B * D * psi * exp
    ]
    eqs.append(gamma * eqs[0] * thinv2 * exp)
    eqs.append(eqs[0] * (-2) * thinv2 * psi)
    eqs.append(du)
    eqs.append(du ** 2)
    return eqs
