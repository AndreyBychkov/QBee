import itertools
from functools import partial, reduce
from random import randrange
from .quadratization import *


def generate_x_plus_1(deg):
    R, x = ring(['x'], QQ)
    return PolynomialSystem([(x + 1) ** deg])


def generate_dense_polynomial(nvars, deg):
    R, *variables = ring([f"x{i}" for i in range(nvars)], QQ)

    def gen_eq_deg(d: int):
        return reduce(add, [randrange(-10, 10) * reduce(mul, p) for p in itertools.product(variables, repeat=d)]) \
               + randrange(-10, 10)

    def gen_eq():
        return reduce(add, [gen_eq_deg(d) for d in range(1, deg + 1)])

    return PolynomialSystem([gen_eq() for _ in range(nvars)])


def generate_circular(deg):
    R, x, y = ring(['x', 'y'], QQ)
    return PolynomialSystem([y ** deg, x ** deg])


def generate_hard(deg):
    R, a, b, c = ring(['a', 'b', 'c'], QQ)
    return PolynomialSystem([c ** deg + a ** 2 * b ** 2 * c ** 3,
                             a ** 2,
                             b ** 2])


def generate_hill(k):
    R, h, i, t = ring(['h', 'i', 't'], QQ)
    return PolynomialSystem([k * i ** 2 * t ** (k - 1),
                             -k * i ** 2 * t ** (k - 1),
                             1])


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
    return PolynomialSystem(system)


def selkov(a, b):
    R, x, y = ring(['x', 'y'], QQ)
    return PolynomialSystem([-x + a * y + x ** 2 * y,
                             b - a * y - x ** 2 * y])
