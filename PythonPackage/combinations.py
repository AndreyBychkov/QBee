import itertools
import sympy as sp

from functools import reduce, partial
from operator import mul
from typing import Tuple, Set, Optional, List


def get_decompositions(monomial: sp.Poly) -> Set[Tuple]:
    r"""
    Returns all decompositions of monomial into submonomials with total degree <= 2.

    Example:
        .. math::
            x^2y^3 \rightarrow (y, xy, xy), (x, xy, y^2), (y, x^2, y^2)

    """
    res = _get_decompositions_rec(monomial, [], set())
    res = map(lambda s: sorted(s, key=sp.polys.orderings.monomial_key('grlex', monomial.free_symbols)), res)
    res = map(lambda b: list(map(lambda m: m.as_expr(), b)), res)  # For convenient representation. Can be removed if necessary.
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
