import sympy as sp

from time import time
from tqdm import tqdm
from combinations import *
from functools import reduce, partial
from typing import List, Tuple
from operator import add
from collections import Counter


def get_possible_replacements(poly_list: List[sp.Expr], excluded_variables: Set[sp.Symbol], count_sorted=True) -> Tuple[sp.Poly]:
    combined_poly = sum(poly_list)
    gens = combined_poly.free_symbols.difference(excluded_variables)
    terms = list(map(lambda m: sp.Poly(m, *gens), sp.Add.make_args(combined_poly.expand())))
    decompositions = list(map(get_all_decompositions, terms))

    possible_replacements = list()
    for dec in decompositions:
        all_replacements = reduce(set.union, map(set, dec))
        possible_replacements += list(filter(lambda p: p.total_degree() > 1, all_replacements))
    possible_replacements = set(possible_replacements)

    if count_sorted:
        decompositions_list = reduce(add, map(list, decompositions))
        decompositions_list = reduce(add, decompositions_list)
        counts = Counter(decompositions_list)
        possible_replacements = sorted(possible_replacements, key=lambda r: counts[r.as_expr()], reverse=True)
    return tuple(possible_replacements)


def polynomial_replace(poly: sp.Expr, old: sp.Expr, new: sp.Expr) -> sp.Expr:
    monomials = sp.Add.make_args(poly)
    return sp.Add(*map(partial(monomial_replace, old=old, new=new), monomials))


def monomial_replace(monomial: sp.Expr, old: sp.Expr, new: sp.Expr) -> sp.Expr:
    quotient, rem = sp.polys.div(monomial, old)
    if rem == 0:
        return quotient * new
    else:
        return monomial


def is_monomial_divisor(numerator: sp.Expr, denominator: sp.Expr) -> bool:
    return sp.gcd(numerator, denominator) == denominator


def sorted_square_first(monomials: List[sp.Poly]) -> List[sp.Poly]:
    return sorted(monomials, key=lambda m: len(m.free_symbols) / m.total_degree())


def reset_progress_bar(pbar: tqdm, value):
    pbar.n = pbar.last_print_n = value
    pbar.start_t = pbar.last_print_t = time()
    pbar.refresh()
