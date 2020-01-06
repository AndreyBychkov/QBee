import sympy as sp

from PythonPackage.AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder
from EquationSystem import EquationSystem
from copy import deepcopy


def is_rational_function(expr: sp.Expr):
    if expr.is_Pow and expr.args1[1] < 0:
        return True
    return False


def preprocess_system(system: EquationSystem) -> None:
    system.expand_equations()


def is_polynomial(system):
    # TODO(can be optimized if we mark polynomized equations)
    non_polynomial_found = None
    for eq in system:
        non_polynomial_found = find_non_polynomial(eq)
    return False if non_polynomial_found else True


def replace_nonlinear_component(system, added_expressions):
    for equation in system:
        non_poly_elem = find_non_polynomial(equation)
        if non_poly_elem:
            new_symbol = symbols_holder.create_symbol()
            added_expressions.append(sp.Eq(new_symbol, non_poly_elem))

            for equation_alt in system:
                equation_alt.subs(non_poly_elem, new_symbol)
            break


def polynomize(eq_system: EquationSystem) -> EquationSystem:
    system = deepcopy(eq_system)
    preprocess_system(system)
    added_expressions = []

    while is_polynomial(system):
        replace_nonlinear_component(system, added_expressions)

    return system
