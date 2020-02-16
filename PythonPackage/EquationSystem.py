import sympy as sp
import random

from functools import reduce
from typing import List, Iterable
from AST_walk import find_non_polynomial
from SymbolsHolder import SymbolsHolder, make_derivative_symbol
from util import polynomial_replace, get_possible_replacements


class EquationSystem:
    def __init__(self, equations: List[sp.Eq],
                 parameter_variables: Iterable[sp.Symbol] = None,
                 input_variables: Iterable[sp.Symbol] = None):
        self._equations = equations.copy()
        self._original_equation_indexes = list(range(len(equations)))
        self._replacement_equations = list()

        self.variables = SymbolsHolder(reduce(set.union, map(lambda e: e.free_symbols, equations)))
        self._parameter_vars = set(parameter_variables) if parameter_variables is not None else set()
        self._input_vars = set(input_variables) if input_variables is not None else set()

    @property
    def equations(self):
        return self._equations

    def replace_expression(self, old: sp.Expr, new: sp.Expr):
        """Replace 'old' expression with 'new' expression for each equation."""
        for i in range(len(self._equations)):
            self._equations[i] = self._equations[i].subs(old, new)

    def replace_subexpression(self, old: sp.Expr, new: sp.Expr):
        """If any expression in system is divisible on 'old', replace it with 'new'"""
        for i, eq in enumerate(self._equations):
            self._equations[i] = sp.Eq(eq.args[0], polynomial_replace(eq.args[1], old, new))

    def expand_equations(self):
        """Apply SymPy 'expand' function to each of equation."""
        for i in range(len(self._equations)):
            self._equations[i] = sp.expand(self._equations[i])

    def is_polynomial(self, mode="original") -> bool:
        """
        Checks if the system is polynomial.

        :param mode: if 'original', checks only original equations of system; if 'full', checks all equations.
        """
        if mode == "original":
            return self._is_polynomial_original()
        elif mode == "full":
            return self._is_polynomial_full()
        else:
            raise ValueError("mode must be 'original' or 'full'.")

    def _is_polynomial_original(self) -> bool:
        for i in self._original_equation_indexes:
            if not self.equations[i].args[1].is_polynomial():
                return False
        return True

    def _is_polynomial_full(self) -> bool:
        for eq in self._equations:
            if not eq.args[1].is_polynomial():
                return False
        return True

    def polynomize(self, mode='algebraic') -> None:
        """
        Transforms the system into polynomial form.

        :param mode: if 'algebraic', adds auxiliary equations in form y = f(x, y);
                     if 'differential', adds auxiliary equations in form y' = f(x, y).
        """
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
        """Calculates Lie derivative using chain rule."""
        result = sp.Integer(0)
        for var in expr.free_symbols.difference(self._parameter_vars).difference(self._input_vars):
            var_diff_eq = list(filter(lambda eq: eq.args[0] == make_derivative_symbol(var), self._equations))[0]
            var_diff = var_diff_eq.args[1]
            result += expr.diff(var) * var_diff
        for input_var in expr.free_symbols.intersection(self._input_vars):
            input_var_dot = make_derivative_symbol(input_var)
            self.variables.add_symbols([input_var_dot])
            result += expr.diff(input_var) * input_var_dot
        return self._replace_from_replacement_equations(result)

    def _replace_from_replacement_equations(self, expr: sp.Expr):
        new_expr = expr.copy()
        for left, right in map(lambda eq: eq.args, self._replacement_equations):
            new_expr = new_expr.subs(right, left)
        return new_expr

    def is_quadratic_linear(self, mode='full') -> bool:
        if mode == 'original':
            return self._is_quadratic_linear_original()
        elif mode == 'full':
            return self._is_quadratic_linear_full()
        else:
            raise ValueError("mode must be 'original' or 'full'.")

    def _is_quadratic_linear_full(self) -> bool:
        for eq in self._equations:
            if not self._is_poly_quadratic_linear(eq.args[1]):
                return False
        return True

    def _is_quadratic_linear_original(self) -> bool:
        for i in self._original_equation_indexes:
            if not self._is_poly_quadratic_linear(self._equations[i].args[1]):
                return False
        return True

    def _is_poly_quadratic_linear(self, poly: sp.Expr) -> bool:
        monomials = sp.Add.make_args(poly)
        for mon in monomials:
            if sp.total_degree(mon) > 2:
                return False
        return True

    def quadratic_linearize(self, mode="algebraic") -> None:
        """
        Transforms the system into quadratic-linear form.

        :param mode: if 'algebraic', adds auxiliary equations in form y = f(x, y);
                     if 'differential', adds auxiliary equations in form y' = f(x, y).
        """
        if not self.is_polynomial():
            raise RuntimeError("System is not polynomized. Polynomize it first.")
        if mode == "algebraic":
            self._quadratic_linearize_algebraic()
        elif mode == "differential":
            self._quadratic_linearize_differential()
        else:
            raise ValueError("mode must be 'algebraic' or 'differential'")

    def _quadratic_linearize_algebraic(self):
        raise NotImplementedError("Algebraic quadratic-linearization is not implemented yet")

    def _quadratic_linearize_differential(self):
        self._quadratic_linearize_diff_rand()

    def _quadratic_linearize_diff_rand(self):
        """Picks replacement in random way."""
        while not self.is_quadratic_linear():
            right_equations = list(map(lambda eq: eq.args[1], self._equations))
            possible_replacements = set(reduce(set.union, map(get_possible_replacements, right_equations)))
            rand_replacement = random.choice(tuple(possible_replacements)).as_expr()

            new_symbol, new_symbol_dot = self.variables.create_symbol_with_derivative()
            self.replace_subexpression(rand_replacement, new_symbol)
            self._replacement_equations.append(sp.Eq(new_symbol, rand_replacement))

            self._equations.append(sp.Eq(new_symbol_dot, self._calculate_Lie_derivative(rand_replacement)).expand())

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))
