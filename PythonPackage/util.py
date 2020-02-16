import sympy as sp

from combinations import get_decompositions
from functools import reduce, partial


def get_possible_replacements(poly: sp.Expr):
    terms = list(map(lambda m: sp.Poly(m), sp.Add.make_args(poly.expand())))
    decompositions = list(map(get_decompositions, terms))

    possible_replacements = list()
    for dec in decompositions:
        all_replacements = reduce(set.union, map(set, dec))
        possible_replacements += list(filter(lambda p: p.total_degree() > 1, map(sp.Poly, all_replacements)))
    return set(possible_replacements)


def polynomial_replace(poly: sp.Expr, old: sp.Expr, new: sp.Expr) -> sp.Expr:
    monomials = sp.Add.make_args(poly)
    return sp.Add(*map(partial(monomial_replace, old=old, new=new), monomials))


def monomial_replace(monomial: sp.Expr, old: sp.Expr, new: sp.Expr) -> sp.Expr:
    quotient, rem = sp.polys.div(monomial, old)
    if rem == 0:
        return quotient * new
    else:
        return monomial
