from sympy.polys.monomials import monomial_deg, monomial_divides
from typing import Tuple, Callable

SelectionStrategy = Callable[['PolynomialSystem'], int]


def empty_score(system) -> int:
    return 1


def default_strategy(system) -> int:
    total_nonsquare = sum([sum(map(abs, m)) for m in system.nonsquares])
    return total_nonsquare + system.dim * len(system.vars)


def aeqd_strategy(system) -> int:
    eq_degs = list(map(lambda mlist: max(map(monomial_deg, mlist)), system.rhs.values()))
    aeqds = list(map(_compute_aeqd, system.nonsquares, [eq_degs, ] * len(system.nonsquares)))
    return sum(aeqds)


def smd_strategy(system) -> int:
    mlist = system.nonsquares
    return sum(map(lambda s: _compute_smd(s, mlist), mlist))


def _compute_aeqd(sub: Tuple[int], eq_degs):
    mon_degs = map(lambda deg: deg + monomial_deg(sub) - 1, eq_degs)
    quad_discrepancies = filter(lambda x: x > 0, map(lambda d: d - 2, mon_degs))
    return sum(quad_discrepancies)


def _compute_smd(sub, mlist: list):
    return (monomial_deg(sub) - 1) * len(list(filter(lambda m: monomial_divides(sub, m), mlist)))
