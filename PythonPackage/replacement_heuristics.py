import sys
import sympy as sp
import random as rand

from structures import *
from functools import partial
from typing import List
from collections import Counter

from util import is_monomial_divisor


def random(system: EquationSystem) -> sp.Expr:
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    rand_replacement = rand.choice(possible_replacements).as_expr()
    return rand_replacement


def count_first(system: EquationSystem) -> sp.Expr:
    possible_replacements = system.get_possible_replacements(True)
    most_frequent_replacement = possible_replacements[0].as_expr()
    return most_frequent_replacement


def sqrt_first(system: EquationSystem) -> sp.Expr:
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    sqrt_replacements = tuple(filter(lambda x: len(x.free_symbols) == 1, possible_replacements))
    if sqrt_replacements:
        return sqrt_replacements[0].as_expr()
    else:
        return possible_replacements[0].as_expr()


def sqrt_count_first(system: EquationSystem) -> sp.Expr:
    possible_replacements = system.get_possible_replacements(count_sorted=True)
    sqrt_replacements = tuple(filter(lambda x: len(x.free_symbols) == 1, possible_replacements))
    if sqrt_replacements:
        return sqrt_replacements[0].as_expr()
    else:
        return possible_replacements[0].as_expr()


def max_replacement_value(system: EquationSystem) -> sp.Expr:
    return max_replacement_value_list(system)[0].as_expr()


def max_replacement_value_list(system: EquationSystem) -> List[sp.Poly]:
    system.update_poly_degrees()  # Optimization can be performed here
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    value_sorted_replacements = sorted(possible_replacements, key=partial(_compute_replacement_value, system))
    return value_sorted_replacements


def _compute_replacement_value(system: EquationSystem, replacement: sp.Poly) -> int:
    used_equations_subs = list(map(make_derivative_symbol, replacement.free_symbols))
    subs_degrees = list(map(lambda subs: system._equations_poly_degrees[subs], used_equations_subs))
    auxiliary_equation_degree = max(subs_degrees)

    return auxiliary_equation_degree


def summary_monomial_degree(system: EquationSystem) -> sp.Expr:
    return summary_monomial_degree(system)[0].as_expr()


def summary_monomial_degree_list(system: EquationSystem) -> List[sp.Poly]:
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    return sorted(possible_replacements, key=partial(_compute_replacement_value_for_all_monomials, system), reverse=False)


def _compute_monomials_affected(system: EquationSystem, replacement: sp.Poly, unique=False) -> int:
    divisible_monomials = filter(partial(is_monomial_divisor, denominator=replacement), system.monomials)
    if unique:
        return len(set(divisible_monomials))
    else:
        return len(list(divisible_monomials))


def _compute_replacement_value_for_all_monomials(system: EquationSystem, replacement: sp.Poly) -> int:
    return replacement.total_degree() * _compute_monomials_affected(system, replacement)


_heuristics_name_to_function = \
    {
        'random'                 : random,
        'sqrt-first'             : sqrt_first,
        'sqrt-count-first'       : sqrt_count_first,
        'replacement-value'      : max_replacement_value,
        'summary-monomial-degree': summary_monomial_degree,
        'default'                : max_replacement_value
    }


def get_heuristics(name: str):
    try:
        return _heuristics_name_to_function[name]
    except KeyError:
        sys.stderr.write("Your heuristics has wrong name. Used default heuristics.")
        return _heuristics_name_to_function['default']
