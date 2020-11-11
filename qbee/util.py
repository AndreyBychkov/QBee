import sympy as sp
import numpy as np
from time import time
from tqdm import tqdm
from typing import Iterable, Collection
from sympy.polys.monomials import monomial_deg
from itertools import combinations
from functools import reduce
from operator import add


def timed(func):
    def wrapper(*args, **kwargs):
        start_time = time()
        res = func(*args, **kwargs)
        end_time = time()
        print()
        print(f"Elapsed time: {np.round(end_time - start_time, 3)}s.")
        return res

    return wrapper


def reset_progress_bar(pbar: tqdm, value):
    pbar.n = pbar.last_print_n = value
    pbar.start_t = pbar.last_print_t = time()
    pbar.refresh()


def refresh_and_close_progress_bars(*pbars: tqdm):
    for bar in pbars:
        bar.refresh()
        bar.close()


def get_decompositions(monomial):
    if len(monomial) == 0:
        return {(tuple(), tuple())}
    result = set()
    prev_result = get_decompositions(tuple(monomial[:-1]))
    for r in prev_result:
        for i in range(monomial[-1] + 1):
            a, b = tuple(list(r[0]) + [i]), tuple(list(r[1]) + [monomial[-1] - i])
            result.add((min(a, b), max(a, b)))
    return result


def symbol_from_derivative(derivative: sp.Symbol) -> sp.Symbol:
    return sp.Symbol(str(derivative).replace(r"\dot ", '', 1))


def poly_to_monomial(poly: sp.Poly) -> sp.Monomial:
    return sp.Monomial(poly.monoms()[0], poly.gens)


def monomial_to_poly(monom: sp.Monomial) -> sp.Poly:
    return sp.Poly(sp.prod([gen ** e for gen, e in zip(monom.gens, monom.exponents)]), monom.gens)


def mlist_to_poly(mlist: Collection[sp.Monomial], gens) -> sp.Poly:
    return sp.Poly(reduce(add, mlist), gens)


def can_substitutions_quadratize(monom: sp.Monomial, subs: Iterable[sp.Monomial]) -> bool:
    gens = tuple(monom.gens)
    var_subs = set(map(lambda var: sp.Monomial(var, gens), gens))
    var_subs.add(sp.Monomial((0,) * len(gens), gens))  # zero degree monomial
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
