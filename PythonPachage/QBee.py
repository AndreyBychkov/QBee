import sympy as sp

from PythonPachage.AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder


def preprocess_system(system: list) -> None:
    for i in range(len(system)):
        system[i] = sp.expand(system[i])


def is_polynomial(system):
    True  # TODO()


def replace_nonlinear_component(system, symbols_holder, added_expressions):
    for equation in system:
        non_poly_elem = find_non_polynomial(equation)
        if non_poly_elem:
            new_symbol = symbols_holder.create_symbol()
            added_expressions.append(sp.Eq(new_symbol, non_poly_elem))

            for equation_alt in system:
                equation_alt.subs(non_poly_elem, new_symbol)
            break


def polynomize(eq_system: list) -> list:
    system = eq_system[:]
    preprocess_system(system)
    symbols_holder = SymbolsHolder(system.free_symbols)
    added_expressions = []

    while is_polynomial(system):
        replace_nonlinear_component(system, symbols_holder, added_expressions)

    return system
