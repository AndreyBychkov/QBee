import sympy as sp

from time import time
from tqdm import tqdm
from functools import partial
from typing import List, Optional, Iterable
from sympy.polys.monomials import monomial_deg
from itertools import combinations


def posify_monomial(monomial: sp.Monomial):
    # TODO: Odd, it was useful with negative substitutions. Has to check this out.
    return monomial if '-' not in str(monomial) else -monomial


def polynomial_subs(poly: sp.Poly, old: sp.Poly, new: sp.Poly) -> sp.Expr:
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


def poly_to_monomial(poly: sp.Poly) -> sp.Monomial:
    return sp.Monomial(poly.monoms()[0], poly.gens)


def monomial_to_poly(monom: sp.Monomial) -> sp.Poly:
    return sp.Poly(sp.prod([gen ** e for gen, e in zip(monom.gens, monom.exponents)]), monom.gens)


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
