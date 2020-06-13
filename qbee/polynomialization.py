import sympy as sp

from copy import deepcopy
from .structures import EquationSystem
from .AST_walk import find_non_polynomial


def polynomialize(system: EquationSystem, mode='differential') -> EquationSystem:
    """
    Transforms the system into polynomial form using variable substitution techniques.

    :param system: non-linear ODEs system
    :param mode: auxiliary equation form.

    Mode
    -----------------
    **algebraic**
        adds auxiliary equations in form y = f(x, y)
    **differential**
         adds auxiliary equations in form y' = f(x, y)

    """
    if mode == 'algebraic':
        return _polynomialize_algebraic(system)
    elif mode == 'differential':
        return _polynomialize_differential(system)
    else:
        raise ValueError("mode must be 'algebraic' or 'differential")


def _polynomialize_algebraic(system: EquationSystem) -> EquationSystem:
    result_system = deepcopy(system)
    while not result_system.is_polynomial():
        _polynomialize_algebraic_iter(result_system)

    return result_system


def _polynomialize_algebraic_iter(system: EquationSystem):
    for eq in system.equations:
        non_poly_elem = find_non_polynomial(eq.args[1])
        if non_poly_elem:
            new_variable = system.variables.create_variable()
            system.algebraic_auxiliary_equation_add(new_variable, non_poly_elem, is_polynomial_substitution=False)
            break


def _polynomialize_differential(system: EquationSystem) -> EquationSystem:
    result_system = deepcopy(system)
    while not result_system.is_polynomial(mode="full"):
        _polynomialize_differential_iter(result_system)

    return result_system


def _polynomialize_differential_iter(system: EquationSystem):
    for eq in system.equations:
        non_poly_elem = find_non_polynomial(eq.args[1])
        if non_poly_elem:

            new_variable = system.variables.create_variable()
            system.differential_auxiliary_equation_add(new_variable, non_poly_elem, is_polynomial_substitution=False)
            break
