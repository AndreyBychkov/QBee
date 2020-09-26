import itertools
import sympy as sp

from sympy.polys.orderings import monomial_key
from functools import reduce, partial
from operator import mul
from typing import Tuple, Set, Optional, List


def get_minimal_decompositions(monomial: sp.Poly) -> Set[Tuple[sp.Poly]]:
    r"""
    Returns decompositions of monomial into submonomials with total degree <= 2.

    Example:
        .. math::
            x^2y^3 \rightarrow (y, xy, xy), (x, xy, y^2), (y, x^2, y^2)

    """
    res = _get_decompositions_rec(monomial, [], set())
    res = res if res is not None else {(monomial,)}
    res = map(lambda s: sorted(s, key=monomial_key('grlex', monomial.free_symbols)), res)
    res = set(tuple(s) for s in res)
    return res


def get_all_decompositions(monomial: sp.Poly) -> Set[Tuple[sp.Poly]]:
    r"""
    Returns all decompositions of monomial into meaning submonomials.

    Example:
        .. math::
            x^2y^3 \rightarrow (x, xy^3), (xy, xy^2), (x^2 y^3), (y, x^2 y^2), (y, x^2, y^2), (x, xy, y^2), (y^2, x^2 y), (y, xy, xy)

    """
    res = _get_complete_decompositions_rec(monomial, [], set())
    res = res if res is not None else {(monomial,)}
    res = set(tuple(s) for s in res)
    return res


def _get_decompositions_rec(monomial: sp.Poly, decomposition: list, result: set) -> Optional[Set]:
    if monomial.total_degree() <= 2:
        return result.add(tuple([monomial] + decomposition))

    for divisor in _get_divisors(monomial):
        _get_decompositions_rec(sp.div(monomial, divisor)[0], decomposition + [divisor], result)
    return result


def _get_divisors(expr: sp.Poly) -> List:
    return _get_possible_squares(expr) + _get_mul_combinations(expr)


def _get_mul_combinations(expr: sp.Poly) -> List:
    res = itertools.combinations(expr.free_symbols, 2)
    res = map(partial(reduce, mul), res)
    res = map(sp.poly, res)
    return list(res)


def _get_possible_squares(expr: sp.Poly) -> List:
    res = filter(lambda var: expr.degree(var) > 1, expr.free_symbols)
    res = map(lambda var: var ** 2, res)
    res = map(sp.poly, res)
    return list(res)


def _get_complete_decompositions_rec(monomial: sp.Poly, decomposition: list, result: set) -> Optional[Set]:
    if decomposition:
        result.add(tuple([monomial] + decomposition))

    if monomial.total_degree() <= 2:
        return

    divisors = list(sp.itermonomials(monomial.gens, monomial.degree_list(), [0] * len(monomial.gens)))
    divisors = list(filter(lambda d: not d.is_Number, divisors))
    divisors = list(map(lambda d: sp.Poly(d, *monomial.gens), divisors))
    divisors = list(filter(lambda d: 1 < d.total_degree() < monomial.total_degree(), divisors))

    for divisor in divisors:
        _get_complete_decompositions_rec(sp.div(monomial, divisor, monomial.gens)[0], decomposition + [divisor], result)

    # Long operation. TODO(try optimize: mb should remove sorting)
    unique_decompositions = set(tuple(sorted(d, key=monomial_key('grlex', monomial.gens))) for d in result)
    return unique_decompositions
