import sys
import sympy as sp
import random as rand

from structures import *
from functools import partial
from typing import List, Tuple
from util import is_monomial_divisor


def random_sorted(system: EquationSystem) -> Tuple[sp.Poly]:
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    rand.shuffle(possible_replacements)
    return possible_replacements


def random(system: EquationSystem) -> sp.Expr:
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    rand_replacement = rand.choice(possible_replacements).as_expr()
    return rand_replacement


def frequent_first_sorted(system: EquationSystem) -> Tuple[sp.Poly]:
    return system.get_possible_replacements(count_sorted=True)


def frequent_first(system: EquationSystem) -> sp.Expr:
    most_frequent_replacement = frequent_first_sorted(system)[0].as_expr()
    return most_frequent_replacement


def free_variables_count_sorted(system: EquationSystem) -> Tuple[sp.Poly]:
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    return tuple(sorted(possible_replacements, key=lambda x: len(x.free_symbols)))


def free_variables_count(system: EquationSystem) -> sp.Expr:
    return free_variables_count_sorted(system)[0].as_expr()


def auxiliary_equation_degree_sorted(system: EquationSystem) -> List[sp.Poly]:
    system.update_poly_degrees()
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    return sorted(possible_replacements, key=partial(_compute_auxiliary_equation_degree, system))


def auxiliary_equation_degree(system: EquationSystem) -> sp.Expr:
    return auxiliary_equation_degree_sorted(system)[0].as_expr()


def _compute_auxiliary_equation_degree(system: EquationSystem, replacement: sp.Poly) -> int:
    used_equations_subs = list(map(make_derivative_symbol, replacement.free_symbols))
    subs_degrees = list(map(lambda subs: system._equations_poly_degrees[subs], used_equations_subs))

    return max(subs_degrees)


def summary_monomial_degree_sorted(system: EquationSystem) -> List[sp.Poly]:
    possible_replacements = system.get_possible_replacements(count_sorted=False)
    return sorted(possible_replacements, key=partial(_compute_replacement_value_for_all_monomials, system), reverse=True)


def summary_monomial_degree(system: EquationSystem) -> sp.Expr:
    return summary_monomial_degree(system)[0].as_expr()


def _compute_monomials_affected(system: EquationSystem, replacement: sp.Poly, unique=False) -> int:
    divisible_monomials = filter(partial(is_monomial_divisor, denominator=replacement), system.monomials)
    if unique:
        return len(set(divisible_monomials))
    else:
        return len(list(divisible_monomials))


def _compute_replacement_value_for_all_monomials(system: EquationSystem, replacement: sp.Poly) -> int:
    return (replacement.total_degree() - 1) * _compute_monomials_affected(system, replacement)


_heuristics_name_to_function = \
    {
        'random'                   : random,
        'frequent-first'           : frequent_first,
        'free-variables-count'     : free_variables_count,
        'auxiliary-equation-degree': auxiliary_equation_degree,
        'summary-monomial-degree'  : summary_monomial_degree,
        'default'                  : summary_monomial_degree
    }

_heuristics_name_to_sorter = \
    {
        'random'                   : random_sorted,
        'frequent-first'           : frequent_first_sorted,
        'free-variables-count'     : free_variables_count_sorted,
        'auxiliary-equation-degree': auxiliary_equation_degree_sorted,
        'summary-monomial-degree'  : summary_monomial_degree_sorted,
        'default'                  : summary_monomial_degree_sorted
    }


def get_heuristics(name: str):
    try:
        return _heuristics_name_to_function[name]
    except KeyError:
        sys.stderr.write("Your heuristics has wrong name. Used default heuristics.")
        return _heuristics_name_to_function['default']


def get_heuristic_sorter(name: str):
    try:
        return _heuristics_name_to_sorter[name]
    except KeyError:
        sys.stderr.write("Your heuristics has wrong name. Used default heuristics.")
        return _heuristics_name_to_sorter['default']
