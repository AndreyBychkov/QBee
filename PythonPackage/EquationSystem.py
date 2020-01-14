import sympy as sp

from functools import reduce
from typing import List
from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder


class EquationSystem:
    def __init__(self, equations: List[sp.Eq]):
        self._equations = equations
        self._original_equation_indexes = list(range(len(equations)))
        self._replacement_equations = list()

        self.variables = SymbolsHolder(reduce(set.union, map(lambda e: e.free_symbols, equations)))

    @property
    def equations(self):
        return self._equations

    def replace_expression(self, old: sp.Expr, new: sp.Expr):
        for i in range(len(self._equations)):
            self._equations[i] = self._equations[i].subs(old, new)

    def expand_equations(self):
        for i in range(len(self._equations)):
            self._equations[i] = sp.expand(self._equations[i])

    def is_polynomial(self):
        for i in self._original_equation_indexes:
            if not self.equations[i].args[1].is_polynomial():
                return False
        return True

    def polynomize(self, mode='algebraic'):
        if mode == 'algebraic':
            self._polynomize_algebraic()
        elif mode == 'differential':
            self._polynomize_differential()
        else:
            raise ValueError("mode must be 'algebraic' or 'differential")

    def _polynomize_algebraic(self):
        while not self.is_polynomial():
            self._replace_algebraic()

    def _replace_algebraic(self):
        for equation in self._equations:
            non_poly_elem = find_non_polynomial(equation)
            if non_poly_elem:
                new_symbol = self.variables.create_symbol()
                self.replace_expression(non_poly_elem, new_symbol)
                self._equations.append(sp.Eq(new_symbol, non_poly_elem))
                break

    def _polynomize_differential(self):
        while not self.is_polynomial():
            self._replace_differential()

    def _replace_differential(self):
        raise NotImplementedError("Differential replace is not implemented.")

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))
