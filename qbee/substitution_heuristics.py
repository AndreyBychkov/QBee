"""
Heuristics choice of monomial substitution for quadratic linearization.
"""

import sys
import sympy as sp
import random as rand

from .structures import *
from functools import partial
from typing import List, Tuple
from .util import is_monomial_divisor


def random_sorted(system: EquationSystem) -> Tuple[sp.Poly]:
    """Shuffles monomial substitutions"""
    possible_substitutions = system.get_possible_substitutions(count_sorted=False)
    return rand.sample(possible_substitutions, len(possible_substitutions))


def random(system: EquationSystem) -> sp.Expr:
    """Picks a random monomial substitution"""
    possible_substitutions = system.get_possible_substitutions(count_sorted=False)
    rand_substitution = rand.choice(possible_substitutions).as_expr()
    return rand_substitution


def frequent_first_sorted(system: EquationSystem) -> Tuple[sp.Poly]:
    """The most frequent monomial substitutions are at the beginning"""
    return system.get_possible_substitutions(count_sorted=True)


def frequent_first(system: EquationSystem) -> sp.Expr:
    """Choose the most frequent monomial substitution"""
    most_frequent_substitution = frequent_first_sorted(system)[0].as_expr()
    return most_frequent_substitution


def free_variables_count_sorted(system: EquationSystem) -> Tuple[sp.Poly]:
    """Monomial substitutions with fewer variables are at the beginning"""
    possible_substitutions = system.get_possible_substitutions(count_sorted=False)
    return tuple(sorted(possible_substitutions, key=lambda x: len(x.free_symbols)))


def free_variables_count(system: EquationSystem) -> sp.Expr:
    """Choose a monomial substitution with the least free variables"""
    return free_variables_count_sorted(system)[0].as_expr()


def auxiliary_equation_degree_sorted(system: EquationSystem) -> List[sp.Poly]:
    """Monomial substitutions with lower generated auxiliary equation degree are at the beginning"""
    system.update_poly_degrees()
    possible_substitutions = system.get_possible_substitutions(count_sorted=False)
    return sorted(possible_substitutions, key=partial(_compute_auxiliary_equation_degree, system))


def auxiliary_equation_degree(system: EquationSystem) -> sp.Expr:
    """Choose a monomial substitution with the least generated auxiliary equation degree"""
    return auxiliary_equation_degree_sorted(system)[0].as_expr()


def _compute_auxiliary_equation_degree(system: EquationSystem, substitution: sp.Poly) -> int:
    used_equations_subs = derivatives(substitution.free_symbols)
    subs_degrees = list(map(lambda subs: system._equations_poly_degrees[subs], used_equations_subs))

    return max(subs_degrees)


def auxiliary_equation_ql_discrepancy_sorted(system: EquationSystem) -> List[sp.Poly]:
    """Monomial substitutions which generated auxiliary equation are closer to quadratic form are at the beginning"""
    system.update_poly_degrees()
    possible_substitutions = system.get_possible_substitutions(count_sorted=False)
    return sorted(possible_substitutions, key=partial(_compute_auxiliary_equation_ql_discrepancy, system))


def auxiliary_equation_ql_discrepancy(system: EquationSystem) -> sp.Expr:
    """Choose a monomial substitution which generated auxiliary equation is closest to quadratic form"""
    return auxiliary_equation_ql_discrepancy_sorted(system)[0].as_expr()


def _compute_auxiliary_equation_ql_discrepancy(system: EquationSystem, substitution: sp.Poly) -> int:
    used_equations_subs = derivatives(substitution.free_symbols)
    monomial_degrees = map(lambda subs: system._equations_poly_degrees[subs] + substitution.total_degree() - 1, used_equations_subs)
    ql_discrepancies = filter(lambda x: x > 0, map(lambda d: d - 2, monomial_degrees))
    return sum(ql_discrepancies)


def summary_monomial_degree_sorted(system: EquationSystem) -> List[sp.Poly]:
    """Monomial substitutions which strongly reduce the degree of the system are at the begging"""
    possible_substitutions = system.get_possible_substitutions(count_sorted=False)
    return sorted(possible_substitutions, key=partial(_compute_substitution_value_for_all_monomials, system), reverse=True)


def summary_monomial_degree(system: EquationSystem) -> sp.Expr:
    """Choose a monomial substitution with maximal reduction of system's degree"""
    return summary_monomial_degree(system)[0].as_expr()


def _compute_monomials_affected(system: EquationSystem, substitution: sp.Poly, unique=False) -> int:
    divisible_monomials = filter(partial(is_monomial_divisor, denominator=substitution), system.monomials)
    if unique:
        return len(set(divisible_monomials))
    else:
        return len(list(divisible_monomials))


def _compute_substitution_value_for_all_monomials(system: EquationSystem, substitution: sp.Poly) -> int:
    return (substitution.total_degree() - 1) * _compute_monomials_affected(system, substitution)


_heuristics_name_to_function = \
    {
        'random'                           : random,
        'frequent-first'                   : frequent_first,
        'free-variables-count'             : free_variables_count,
        'auxiliary-equation-degree'        : auxiliary_equation_degree,
        'auxiliary-equation-ql-discrepancy': auxiliary_equation_ql_discrepancy,
        'summary-monomial-degree'          : summary_monomial_degree,
        'default'                          : summary_monomial_degree
    }

_heuristics_name_to_sorter = \
    {
        'random'                           : random_sorted,
        'frequent-first'                   : frequent_first_sorted,
        'free-variables-count'             : free_variables_count_sorted,
        'auxiliary-equation-degree'        : auxiliary_equation_degree_sorted,
        'auxiliary-equation-ql-discrepancy': auxiliary_equation_ql_discrepancy_sorted,
        'summary-monomial-degree'          : summary_monomial_degree_sorted,
        'default'                          : summary_monomial_degree_sorted
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
