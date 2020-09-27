import sympy as sp

from sympy.polys.monomials import monomial_deg
from sympy.polys.orderings import monomial_key
from typing import List, Iterable, Collection
from operator import add
from functools import reduce


def get_substitutions(poly_list: List[sp.Poly], gens: List[sp.Symbol] = None) -> Collection[sp.Monomial]:
    if gens is None:
        gens = poly_list[0].gens

    gens_greatest_monoms = _get_greatest_monomials(poly_list, gens)
    gcd_monom = gcd(*gens_greatest_monoms)
    subs = sorted(sp.itermonomials(gens, gcd_monom.exponents, [0] * len(gens)), key=monomial_key('grlex', gens)) # slice in order to remove monomials of degree 0 and 1
    for gr_monom in gens_greatest_monoms:
        gr_monom_subs = list(sp.itermonomials(gens, gr_monom.exponents, gcd_monom.exponents))[:-1]
        subs += gr_monom_subs
    subs = [sp.Monomial(s, gens) for s in set(subs)]
    return list(filter(lambda m: monomial_deg(m) > 1, subs))
    # Todo: Last monomials by grlex still remain. It could affect quadratization algorithm. Need to fix this if overhead is not so big.


def _get_greatest_monomials(poly_list: List[sp.Poly], gens: List[sp.Symbol]) -> Iterable[sp.Monomial]:
    monoms = list(set(reduce(add, map(sp.Poly.monoms, poly_list))))
    greatest_monomials_lists = [maxes(monoms, key=lambda a: a[i]) for i, _ in enumerate(gens)]
    greatest_monomials = [max(mon_lst, key=monomial_deg) for mon_lst in greatest_monomials_lists]
    return [sp.Monomial(m, gens) for m in greatest_monomials]


def gcd(*args: sp.Monomial) -> sp.Monomial:
    return sp.Monomial([min(m_i) for m_i in zip(*args)])


def maxes(a, key=None):
    if key is None:
        key = lambda x: x
    m, max_list = key(a[0]), []
    for s in a:
        k = key(s)
        if k > m:
            m, max_list = k, [s]
        elif k == m:
            max_list.append(s)
    return max_list
