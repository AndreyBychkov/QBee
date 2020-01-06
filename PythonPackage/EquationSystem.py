import sympy as sp

from functools import reduce
from typing import List
from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder


class EquationSystem:
    def __init__(self, equations: List[sp.Eq]):
        self._equations = equations
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
        non_polynomial_found = None
        for eq in self._equations:
            non_polynomial_found = find_non_polynomial(eq)
        return False if non_polynomial_found else True

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))
