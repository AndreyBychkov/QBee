import sympy as sp

from time import time
from tqdm import tqdm
from .combinations import *
from itertools import combinations
from functools import reduce, partial
from typing import List, Tuple, Iterable
from collections import Counter


def get_possible_substitutions(poly_list: List[sp.Expr], gens: Set[sp.Symbol], count_sorted=True) -> Tuple[sp.Poly]:
    decompositions = get_monomial_decompositions(poly_list, gens)
    return get_possible_substitutions_from_decompositions(decompositions, count_sorted)


def get_monomial_decompositions(poly_list: List[sp.Expr], gens: Set[sp.Symbol]) -> List[Set[Tuple[sp.Poly]]]:
    polys = list(map(lambda m: sp.Poly(m, *gens).expand(), poly_list))
    terms = reduce(lambda a, b: a + b, map(lambda p: sp.Add.make_args(p), polys))
    return list(map(get_all_decompositions, terms))


def unify_poly_field(poly_list: List[sp.Poly]) -> List[sp.Poly]:
    gens = set(reduce(lambda a, b: a + b, map(lambda p: p.gens, poly_list)))
    return list(map(lambda p: sp.Poly(p, *gens), poly_list))


def get_possible_substitutions_from_decompositions(decompositions: List[Set[Tuple[sp.Poly]]], count_sorted=False) -> Tuple[sp.Poly]:
    possible_substitutions = list()
    for dec in decompositions:
        all_substitutions = reduce(set.union, map(set, dec))
        possible_substitutions += list(filter(lambda p: p.total_degree() > 1, all_substitutions))
    possible_substitutions = set(map(posify_monomial, possible_substitutions))

    if count_sorted:
        decompositions_list = sum(map(list, decompositions))
        decompositions_list = sum(decompositions_list)
        counts = Counter(decompositions_list)
        possible_substitutions = sorted(possible_substitutions, key=lambda r: counts[r.as_expr()], reverse=True)
    return tuple(possible_substitutions)


def posify_monomial(monomial: sp.Monomial):
    return monomial if '-' not in str(monomial) else -monomial


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


def reset_progress_bar(pbar: tqdm, value):
    pbar.n = pbar.last_print_n = value
    pbar.start_t = pbar.last_print_t = time()
    pbar.refresh()


def refresh_and_close_progress_bars(*pbars: tqdm):
    for bar in pbars:
        bar.refresh()
        bar.close()


def symbol_from_derivative(derivative: sp.Symbol) -> sp.Symbol:
    return sp.Symbol(str(derivative).replace(r"\dot ", '', 1))


def can_substitutions_quadratize(monom: sp.Monomial, subs: Iterable[sp.Monomial]) -> bool:
    gens = monom.gens
    expanded_subs = set(map(lambda var: sp.Monomial(var, monom.gens), gens)).union(subs)
    return _can_quad_2(monom, expanded_subs) or _can_quad_1(monom, expanded_subs)


def _can_quad_2(monom: sp.Monomial, subs: Iterable[sp.Monomial]) -> bool:
    for sub in subs:
        if monom == sub ** 2:
            return True
    return False


def _can_quad_1(monom: sp.Monomial, subs: Iterable[sp.Monomial]) -> bool:
    for left, right in combinations(subs, 2):
        if monom == left * right:
            return True
    return False
