from sympy.polys.monomials import monomial_deg, monomial_divides
from typing import Tuple, Callable, Iterable, Any

from .util import *

GenerationStrategy = Callable[['PolynomialSystem'], Iterable[Any]]
Scoring = Callable[['PolynomialSystem'], int]

# Standard generation strategy with different scoring functions

def default_generation(system):
    if len(system.nonsquares) == 0:
        return list()
    return get_decompositions(system.get_smallest_nonsquare())

# Different scoring functions

def empty_score(system) -> int:
    return 1

def default_scoring(system) -> int:
    total_nonsquare = sum([sum(map(abs, m)) for m in system.nonsquares])
    return total_nonsquare + system.dim * len(system.vars)


def aeqd_scoring(system) -> int:
    eq_degs = list(map(lambda mlist: max(map(monomial_deg, mlist)), system.rhs.values()))
    aeqds = list(map(_compute_aeqd, system.nonsquares, [eq_degs, ] * len(system.nonsquares)))
    return sum(aeqds)


def smd_scoring(system) -> int:
    mlist = system.nonsquares
    return sum(map(lambda s: _compute_smd(s, mlist), mlist))


def _compute_scoring(sub: Tuple[int], eq_degs):
    mon_degs = map(lambda deg: deg + monomial_deg(sub) - 1, eq_degs)
    quad_discrepancies = filter(lambda x: x > 0, map(lambda d: d - 2, mon_degs))
    return sum(quad_discrepancies)


def _compute_scoring(sub, mlist: list):
    return (monomial_deg(sub) - 1) * len(list(filter(lambda m: monomial_divides(sub, m), mlist)))
