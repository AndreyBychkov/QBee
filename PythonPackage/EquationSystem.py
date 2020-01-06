import sympy as sp

from functools import reduce
from typing import List
from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder


class EquationSystem:
    def __init__(self, equations: List[sp.Eq]):
        self.original_equations = equations
        self.added_equations = list()
        self.replacement_equations = list()

        self.variables = SymbolsHolder(reduce(set.union, map(lambda e: e.free_symbols, equations)))

    def get_equations(self):
        return self.original_equations + self.added_equations

    def replace_expression(self, old: sp.Expr, new: sp.Expr):
        for i in range(len(self.original_equations)):
            self.original_equations[i] = self.original_equations[i].subs(old, new)

    def expand_equations(self):
        for i in range(len(self.original_equations)):
            self.original_equations[i] = sp.expand(self.original_equations[i])

        for i in range(len(self.added_equations)):
            self.added_equations[i] = sp.expand(self.added_equations[i])

    def is_polynomial(self):
        non_polynomial_found = None
        for eq in self.original_equations + self.added_equations:
            non_polynomial_found = find_non_polynomial(eq)
        return False if non_polynomial_found else True

    def __len__(self):
        return self.original_equations.__len__() + self.added_equations.__len__()

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self.get_equations()))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self.get_equations()))
