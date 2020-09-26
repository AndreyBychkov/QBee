"""
Heuristics choice of monomial substitution for quadratization.
"""

import sys
import sympy as sp
import random as rand

from .structures import PolynomialSystem, derivatives
from functools import partial
from typing import List, Tuple, Sequence, Callable
from .util import is_monomial_divisor


def none_sorted(system: PolynomialSystem) -> Sequence[sp.Poly]:
    """Does not perform any sorting"""
    return system.get_possible_substitutions()


def none_subs_sorted(system: PolynomialSystem, substitutions: Sequence[sp.Poly]) -> Sequence[sp.Poly]:
    return substitutions


def choose_first(system: PolynomialSystem) -> sp.Poly:
    return none_sorted(system)[0].as_expr()


def random_sorted(system: PolynomialSystem) -> Sequence[sp.Poly]:
    """Shuffles monomial substitutions"""
    possible_substitutions = system.get_possible_substitutions()
    return random_subs_sorted(system, possible_substitutions)


def random_subs_sorted(system: PolynomialSystem, substitutions: Sequence[sp.Poly]) -> Sequence[sp.Poly]:
    return rand.sample(substitutions, len(substitutions))


def random(system: PolynomialSystem) -> sp.Poly:
    """Picks a random monomial substitution"""
    possible_substitutions = system.get_possible_substitutions()
    rand_substitution = rand.choice(possible_substitutions).as_expr()
    return rand_substitution


def free_variables_count_sorted(system: PolynomialSystem) -> Sequence[sp.Poly]:
    """Monomial substitutions with fewer variables are at the beginning"""
    possible_substitutions = system.get_possible_substitutions()
    return free_variables_count_subs_sorted(system, possible_substitutions)


def free_variables_count_subs_sorted(system: PolynomialSystem, substitutions: Sequence[sp.Poly]) -> Sequence[sp.Poly]:
    """Monomial substitutions with fewer variables are at the beginning"""
    return tuple(sorted(substitutions, key=lambda x: len(x.free_symbols)))


def free_variables_count(system: PolynomialSystem) -> sp.Poly:
    """Choose a monomial substitution with the least free variables"""
    return free_variables_count_sorted(system)[0].as_expr()


def auxiliary_equation_degree_sorted(system: PolynomialSystem) -> Sequence[sp.Poly]:
    """Monomial substitutions with lower generated auxiliary equation degree are at the beginning"""
    possible_substitutions = system.get_possible_substitutions()
    return auxiliary_equation_degree_subs_sorted(system, possible_substitutions)


def auxiliary_equation_degree_subs_sorted(system: PolynomialSystem, substitutions: Sequence[sp.Poly]) -> Sequence[sp.Poly]:
    """Monomial substitutions with lower generated auxiliary equation degree are at the beginning"""
    system.update_poly_degrees()
    return sorted(substitutions, key=partial(_compute_auxiliary_equation_degree, system))


def auxiliary_equation_degree(system: PolynomialSystem) -> sp.Poly:
    """Choose a monomial substitution with the least generated auxiliary equation degree"""
    return auxiliary_equation_degree_sorted(system)[0].as_expr()


def _compute_auxiliary_equation_degree(system: PolynomialSystem, substitution: sp.Poly) -> int:
    used_equations_subs = derivatives(substitution.free_symbols)
    subs_degrees = list(map(lambda subs: system._equations_poly_degrees[subs], used_equations_subs))
    return max(subs_degrees)


def auxiliary_equation_quadratic_discrepancy_sorted(system: PolynomialSystem) -> Sequence[sp.Poly]:
    """Monomial substitutions which generated auxiliary equation are closer to quadratic form are at the beginning"""
    possible_substitutions = system.get_possible_substitutions()
    return auxiliary_equation_quadratic_discrepancy_subs_sorted(system, possible_substitutions)


def auxiliary_equation_quadratic_discrepancy_subs_sorted(system: PolynomialSystem, substitutions: Sequence[sp.Poly]) -> Sequence[sp.Poly]:
    """Monomial substitutions which generated auxiliary equation are closer to quadratic form are at the beginning"""
    system.update_poly_degrees()
    return sorted(substitutions, key=partial(_compute_auxiliary_equation_quadratic_discrepancy, system))


def auxiliary_equation_quadratic_discrepancy(system: PolynomialSystem) -> sp.Poly:
    """Choose a monomial substitution which generated auxiliary equation is closest to quadratic form"""
    return auxiliary_equation_quadratic_discrepancy_sorted(system)[0].as_expr()


def _compute_auxiliary_equation_quadratic_discrepancy(system: PolynomialSystem, substitution: sp.Poly) -> int:
    used_equations_subs = derivatives(substitution.free_symbols.intersection(system.variables.free))
    monomial_degrees = map(lambda subs: system._equations_poly_degrees[subs] + substitution.total_degree() - 1, used_equations_subs)
    quad_discrepancies = filter(lambda x: x > 0, map(lambda d: d - 2, monomial_degrees))
    return sum(quad_discrepancies)


def summary_monomial_degree_sorted(system: PolynomialSystem) -> Sequence[sp.Poly]:
    """Monomial substitutions which strongly reduce the degree of the system are at the begging"""
    possible_substitutions = system.get_possible_substitutions()
    return summary_monomial_degree_subs_sorted(system, possible_substitutions)


def summary_monomial_degree_subs_sorted(system: PolynomialSystem, substitutions: Sequence[sp.Poly]) -> Sequence[sp.Poly]:
    """Monomial substitutions which strongly reduce the degree of the system are at the begging"""
    return sorted(substitutions, key=partial(_compute_substitution_value_for_all_monomials, system), reverse=True)


def summary_monomial_degree(system: PolynomialSystem) -> sp.Poly:
    """Choose a monomial substitution with maximal reduction of system's degree"""
    return summary_monomial_degree(system)[0].as_expr()


def _compute_monomials_affected(system: PolynomialSystem, substitution: sp.Poly, unique=False) -> int:
    divisible_monomials = filter(partial(is_monomial_divisor, denominator=substitution), system.monomials)
    if unique:
        return len(set(divisible_monomials))
    else:
        return len(list(divisible_monomials))


def _compute_substitution_value_for_all_monomials(system: PolynomialSystem, substitution: sp.Poly) -> int:
    return (substitution.total_degree() - 1) * _compute_monomials_affected(system, substitution)


_heuristics_name_to_function = \
    {
        'none'                                    : choose_first,
        'random'                                  : random,
        'free-variables-count'                    : free_variables_count,
        'auxiliary-equation-degree'               : auxiliary_equation_degree,
        'auxiliary-equation-quadratic-discrepancy': auxiliary_equation_quadratic_discrepancy,
        'summary-monomial-degree'                 : summary_monomial_degree,
        'FVC'                                     : free_variables_count,
        'AED'                                     : auxiliary_equation_degree,
        'AEQD'                                    : auxiliary_equation_quadratic_discrepancy,
        'SMD'                                     : summary_monomial_degree,
        'default'                                 : summary_monomial_degree
    }

_heuristics_name_to_sorter = \
    {
        'none'                                    : none_sorted,
        'random'                                  : random_sorted,
        'free-variables-count'                    : free_variables_count_sorted,
        'auxiliary-equation-degree'               : auxiliary_equation_degree_sorted,
        'auxiliary-equation-quadratic-discrepancy': auxiliary_equation_quadratic_discrepancy_sorted,
        'summary-monomial-degree'                 : summary_monomial_degree_sorted,
        'FVC'                                     : free_variables_count_sorted,
        'AED'                                     : auxiliary_equation_degree_sorted,
        'AEQD'                                    : auxiliary_equation_quadratic_discrepancy_sorted,
        'SMD'                                     : summary_monomial_degree_sorted,
        'default'                                 : summary_monomial_degree_sorted
    }

_heuristics_name_to_subs_sorter = \
    {
        'none'                                    : none_subs_sorted,
        'random'                                  : random_subs_sorted,
        'free-variables-count'                    : free_variables_count_subs_sorted,
        'auxiliary-equation-degree'               : auxiliary_equation_degree_subs_sorted,
        'auxiliary-equation-quadratic-discrepancy': auxiliary_equation_quadratic_discrepancy_subs_sorted,
        'summary-monomial-degree'                 : summary_monomial_degree_subs_sorted,
        'FVC'                                     : free_variables_count_subs_sorted,
        'AED'                                     : auxiliary_equation_degree_subs_sorted,
        'AEQD'                                    : auxiliary_equation_quadratic_discrepancy_subs_sorted,
        'SMD'                                     : summary_monomial_degree_subs_sorted,
        'default'                                 : summary_monomial_degree_subs_sorted
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


def get_heuristics_substitution_sorted(name: str):
    try:
        return _heuristics_name_to_subs_sorter[name]
    except KeyError:
        sys.stderr.write("Your heuristics has wrong name. Used default heuristics.")
        return _heuristics_name_to_sorter['default']
