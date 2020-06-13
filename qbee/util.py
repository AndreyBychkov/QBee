import sympy as sp

from time import time
from tqdm import tqdm
from .combinations import *
from functools import reduce, partial
from typing import List, Tuple
from operator import add
from collections import Counter


def get_possible_substitutions(poly_list: List[sp.Expr], gens: Set[sp.Symbol], count_sorted=True) -> Tuple[sp.Poly]:
    combined_poly = sum(poly_list)
    terms = list(map(lambda m: sp.Poly(m, *gens), sp.Add.make_args(combined_poly.expand())))
    decompositions = list(map(get_all_decompositions, terms))

    possible_substitutions = list()
    for dec in decompositions:
        all_substitutions = reduce(set.union, map(set, dec))
        possible_substitutions += list(filter(lambda p: p.total_degree() > 1, all_substitutions))
    possible_substitutions = set(possible_substitutions)

    if count_sorted:
        decompositions_list = reduce(add, map(list, decompositions))
        decompositions_list = reduce(add, decompositions_list)
        counts = Counter(decompositions_list)
        possible_substitutions = sorted(possible_substitutions, key=lambda r: counts[r.as_expr()], reverse=True)
    return tuple(possible_substitutions)


def polynomial_subs(poly: sp.Expr, old: sp.Expr, new: sp.Expr) -> sp.Expr:
    monomials = sp.Add.make_args(poly)
    return sp.Add(*map(partial(monomial_subs, old=old, new=new), monomials))


def monomial_subs(monomial: sp.Expr, old: sp.Expr, new: sp.Expr) -> sp.Expr:
    quotient, rem = sp.polys.div(monomial, old)
    if rem == 0:
        return quotient * new
    else:
        return monomial


def is_monomial_divisor(numerator: sp.Expr, denominator: sp.Expr) -> bool:
    return sp.gcd(numerator.as_expr(), denominator.as_expr()) == denominator.as_expr()


def sorted_square_first(monomials: List[sp.Poly]) -> List[sp.Poly]:
    return sorted(monomials, key=lambda m: len(m.free_symbols) / m.total_degree())


def reset_progress_bar(pbar: tqdm, value):
    pbar.n = pbar.last_print_n = value
    pbar.start_t = pbar.last_print_t = time()
    pbar.refresh()


def refresh_and_close_progress_bars(*pbars: tqdm):
    for bar in pbars:
        bar.refresh()
        bar.close()
