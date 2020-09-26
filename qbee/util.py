import sympy as sp

from time import time
from tqdm import tqdm
from .combinations import *
from sympy.polys.monomials import monomial_deg
from itertools import combinations
from functools import reduce, partial
from typing import List, Tuple, Iterable, Collection
from collections import Counter


def get_possible_substitutions(poly_list: List[sp.Expr], gens: Set[sp.Symbol], count_sorted=True) -> Tuple[sp.Poly]:
    decompositions = get_monomial_decompositions(poly_list, gens)
    return get_possible_substitutions_from_decompositions(decompositions, count_sorted)


def get_monomial_decompositions(poly_list: List[sp.Expr], gens: Set[sp.Symbol]) -> List[Set[Tuple[sp.Poly]]]:
    terms = reduce(lambda a, b: a + b, map(lambda p: sp.Add.make_args(p.as_expr()), poly_list))
    return list(map(get_all_decompositions, [sp.Poly(t, *gens) for t in terms]))


def unify_poly_field(poly_list: List[sp.Poly], gens: Collection[sp.Symbol] = None) -> List[sp.Poly]:
    if gens is None:
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
    gens = tuple(monom.gens)
    var_subs = set(map(lambda var: sp.Monomial(var, gens), gens))
    var_subs.add(sp.Monomial((0,) * len(gens), gens))  # zero degree monomial
    # unified_gens_subs = set(map(lambda sub: unify_monom_to_gens(sub, gens), subs))
    expanded_subs = var_subs.union(subs)
    return _can_quad_0(monom) or _can_quad_2(monom, expanded_subs) or _can_quad_1(monom, expanded_subs)


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


def _can_quad_0(monom: sp.Monomial) -> bool:
    return monomial_deg(monom) == 0


def poly_to_monomial(poly: sp.Poly) -> sp.Monomial:
    return sp.Monomial(poly.monoms()[0], poly.gens)


def unify_monom_to_gens(in_monom: sp.Monomial, gens: List[sp.Symbol]):
    res = list()
    for gen in gens:
        gen_index = index_or_None(in_monom.gens, gen)
        if gen_index is not None:
            res.append(in_monom.exponents[gen_index])
        else:
            res.append(0)
    return sp.Monomial(res, gens)


def index_or_None(it: Iterable, elem) -> Optional[int]:
    for i, it_elem in enumerate(it):
        if it_elem == elem:
            return i
    return None


def gcd(*args: sp.Monomial) -> sp.Monomial:
    return sp.Monomial([min(m_i) for m_i in zip(*args)])
