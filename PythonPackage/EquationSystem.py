import sympy as sp

from functools import reduce
from typing import List
from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder, make_derivative_symbol


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

    def is_polynomial(self, mode="original"):
        if mode == "original":
            return self._is_polynomial_original()
        elif mode == "full":
            return self._is_polynomial_full()
        else:
            raise ValueError("mode must be 'original' or 'full'.")

    def _is_polynomial_original(self):
        for i in self._original_equation_indexes:
            if not self.equations[i].args[1].is_polynomial():
                return False
        return True

    def _is_polynomial_full(self):
        for eq in self._equations:
            if not eq.args[1].is_polynomial():
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
        for i in self._original_equation_indexes:
            non_poly_elem = find_non_polynomial(self.equations[i].args[1])
            if non_poly_elem:
                new_symbol = self.variables.create_symbol()
                self.replace_expression(non_poly_elem, new_symbol)
                self._replacement_equations.append(sp.Eq(new_symbol, non_poly_elem))
                self._equations.append(sp.Eq(new_symbol, non_poly_elem))
                break

    def _polynomize_differential(self):
        while not self.is_polynomial(mode="full"):
            self._replace_differential()

    def _replace_differential(self):
        for eq in self._equations:
            non_poly_elem = find_non_polynomial(eq.args[1])
            if non_poly_elem:
                new_symbol, new_symbol_der = self.variables.create_symbol_with_derivative()
                self.replace_expression(non_poly_elem, new_symbol)
                self._replacement_equations.append(sp.Eq(new_symbol, non_poly_elem))
                self._equations.append(sp.Eq(new_symbol_der, self._calculate_Lie_derivative(non_poly_elem)))
                break

    def _calculate_Lie_derivative(self, expr: sp.Expr):
        der = expr.diff(*expr.free_symbols)
        for var in expr.free_symbols:
            var_diff_eq = list(filter(lambda eq: eq.args[0] == make_derivative_symbol(var), self._equations))[0]
            var_diff = var_diff_eq.args[1]
            der *= var_diff
        return self._replace_from_list(der)

    def _replace_from_list(self, expr: sp.Expr):
        new_expr = expr.copy()
        for left, right in map(lambda eq: eq.args, self._replacement_equations):
            new_expr = new_expr.subs(right, left)
        return new_expr

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))
