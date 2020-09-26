import sympy as sp

from sympy.polys.orderings import monomial_key
from typing import Set, List, Iterable, Collection
from operator import add
from functools import reduce


def get_substitutions(poly_list: List[sp.Poly], gens: List[sp.Symbol] = None) -> Collection[sp.Monomial]:
    if gens is None:
        gens = poly_list[0].gens

    gens_greatest_monoms = _get_greatest_monomials(poly_list, gens)
    gcd_monom = gcd(*gens_greatest_monoms)
    subs = sorted(sp.itermonomials(gens, gcd_monom.exponents, [0] * len(gens)), key=monomial_key('grlex', gens))[
           len(gens) + 1:]  # slice in order to remove monomials of degree 0 and 1
    for gr_monom in gens_greatest_monoms:
        gr_monom_subs = list(sp.itermonomials(gens, gr_monom.exponents, gcd_monom.exponents))[:-1]
        subs += gr_monom_subs
    return [sp.Monomial(s, gens) for s in set(subs)]
    # Todo: Last monomials by grlex still remain. It could affect quadratization algorithm. Need to fix this if overhead is not so big.


def _get_greatest_monomials(poly_list: List[sp.Poly], gens: List[sp.Symbol]) -> Iterable[sp.Monomial]:
    monoms = set(reduce(add, map(sp.Poly.monoms, poly_list)))
    greatest_monomials = [max(monoms, key=lambda a: a[i]) for i, _ in enumerate(gens)]
    return [sp.Monomial(m, gens) for m in greatest_monomials]


def gcd(*args: sp.Monomial) -> sp.Monomial:
    return sp.Monomial([min(m_i) for m_i in zip(*args)])
