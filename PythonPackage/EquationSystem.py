import os
import random
import hashlib
import sympy as sp
import pandas as pd

from tqdm import tqdm
from copy import copy, deepcopy
from queue import Queue
from functools import reduce
from typing import List, Iterable, Optional, Callable
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

    @property
    def equations_hash(self):
        return hashlib.md5(str(self._equations).encode('utf-8')).digest()

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

    def polynomialize(self, mode='differential') -> None:
        """
        Transforms the system into polynomial form using variable replacement techniques.

        :param mode: auxiliary equation form.

        Mode
        -----------------
        **algebraic**
            adds auxiliary equations in form y = f(x, y)
        **differential**
             adds auxiliary equations in form y' = f(x, y)

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

    def quadratic_linearized(self, mode="heuristic", auxiliary_eq_type="differential", heuristics='sqrt-count-first', debug=None, log_file=None):
        """
        Transforms the system into quadratic-linear form using variable replacement technique.

        :param mode: use 'optimal' to find optimal transformation.
        :param auxiliary_eq_type: auxiliary equation form.
        :param heuristics: next replacement choice method.
        :param debug: printing mode while quadratic linearization is performed.
        :param log_file: output file for evaluation logging. Must be in 'csv' format.
        :returns: quadtaric-linearizaed system
        :rtype: EquationSystem

        Mode
        -----------------
        **optimal**
            find optimal transformation. The most time-consuming mode;
        **heuristic**
            find sub-optimal transformation. Works much faster than 'optimal'. You can choose heuristics in 'heuristics' parameter;

        Auxiliary equations type
        -----------------
        **algebraic**
            adds auxiliary equations in form y = f(x, y)
        **differential**
             adds auxiliary equations in form y' = f(x, y)

        Heuristics
        -----------------
        **random**
            choose next possible replacement in random way;
        **count-first**
            choose most frequent possible replacement as the next one;
        **sqrt-first**
            choose next possible replacement within variables' squares in random way;
        **sqrt-count-first**
             choose most frequent square replacement as the next one;

        Debug
        ---------------
        **None** or **silent**
            prints nothing;
        **info**
            prints replacement for each iteration;
        **debug**
            prints equations in system with replacement for each iteration;

        """
        if not self.is_polynomial():
            raise RuntimeError("System is not polynomialized. Polynomize it first.")
        if mode == 'optimal':
            if auxiliary_eq_type == 'algebraic':
                return self._ql_algebraic_optimal()
            elif auxiliary_eq_type == 'differential':
                return self._ql_differential_optimal()
            else:
                raise ValueError("auxiliary_eq_type must be 'algebraic' or 'differential'")
        elif mode == 'heuristic':
            if auxiliary_eq_type == "algebraic":
                return self._quadratic_linearize_algebraic(debug=debug)
            elif auxiliary_eq_type == "differential":
                return self._quadratic_linearize_differential(heuristics, debug=debug, log_file=log_file)
            else:
                raise ValueError("auxiliary_eq_type must be 'algebraic' or 'differential'")
        else:
            raise ValueError("mode must be 'optimal' or 'heuristic'")

    def _quadratic_linearize_algebraic(self, debug=None):
        raise NotImplementedError("Algebraic quadratic-linearization is not implemented yet")

    def _quadratic_linearize_differential(self, method: str, debug: Optional[str] = None, log_file: Optional[str] = None):
        log_rows_list = list()
        new_system = deepcopy(self)
        while not new_system.is_quadratic_linear():
            hash_before, hash_after, replacement = new_system._ql_differential_iter(method)

            new_system._debug_system_print(debug)
            if log_file:
                new_system._ql_log_append(log_rows_list, hash_before, hash_after, replacement)

        if log_file:
            log_df = pd.DataFrame(log_rows_list)
            log_df.to_csv(log_file, index=False)

        if not (debug is None or debug == 'silent'):
            print('-' * 100)

        return new_system

    def _ql_differential_iter(self, method: str):
        hash_before = self.equations_hash

        replacement = self._ql_choice_method_name_to_function(method)()
        new_symbol, new_symbol_dot = self.variables.create_symbol_with_derivative()
        self.replace_subexpression(replacement, new_symbol)
        self._replacement_equations.append(sp.Eq(new_symbol, replacement))

        self._equations.append(sp.Eq(new_symbol_dot, self._calculate_Lie_derivative(replacement)).expand())
        hash_after = self.equations_hash

        return hash_before, hash_after, replacement

    def _ql_log_append(self, row_list: List, hash_before, hash_after, replacement):
        row_list.append({'from': hash_before, 'name': hash_after, 'replacement': replacement})

    def _ql_choice_method_name_to_function(self, method):
        if method == 'random':
            return self._ql_random_choice
        elif method == 'count-first':
            return self._ql_count_first_choice
        elif method == 'sqrt-first':
            return self._ql_sqrt_first_choice
        elif method == 'sqrt-count-first':
            return self._ql_sqrt_count_first_choice
        else:
            raise ValueError("Replacement method has wrong name.")

    def _ql_random_choice(self):
        right_equations = list(map(lambda eq: eq.args[1], self._equations))
        possible_replacements = get_possible_replacements(right_equations, count_sorted=False)
        rand_replacement = random.choice(possible_replacements).as_expr()
        return rand_replacement

    def _ql_count_first_choice(self):
        right_equations = list(map(lambda eq: eq.args[1], self._equations))
        possible_replacements = get_possible_replacements(right_equations, count_sorted=True)
        most_frequent_replacement = possible_replacements[0].as_expr()
        return most_frequent_replacement

    def _ql_sqrt_first_choice(self):
        right_equations = list(map(lambda eq: eq.args[1], self._equations))
        possible_replacements = get_possible_replacements(right_equations, count_sorted=False)
        sqrt_replacements = tuple(filter(lambda x: len(x.free_symbols) == 1, possible_replacements))
        if sqrt_replacements:
            return sqrt_replacements[0].as_expr()
        else:
            return possible_replacements[0].as_expr()

    def _ql_sqrt_count_first_choice(self):
        right_equations = list(map(lambda eq: eq.args[1], self._equations))
        possible_replacements = get_possible_replacements(right_equations, count_sorted=True)
        sqrt_replacements = tuple(filter(lambda x: len(x.free_symbols) == 1, possible_replacements))
        if sqrt_replacements:
            return sqrt_replacements[0].as_expr()
        else:
            return possible_replacements[0].as_expr()

    def _ql_algebraic_optimal(self):
        raise NotImplementedError("Algebraic optimal quadratic-linearization is not implemented yet. ")

    def _ql_differential_optimal(self):
        system_queue = Queue()
        system_queue.put(self, block=True)
        ql_reached = False
        counter = 0
        while not ql_reached:
            if counter % 50 == 0:
                os.system('cls')
                print(f"Systems processed: {counter}")
            curr_system = system_queue.get()
            if curr_system.is_quadratic_linear():
                return curr_system
            right_equations = list(map(lambda eq: eq.args[1], curr_system.equations))
            possible_replacements = get_possible_replacements(right_equations, count_sorted=False)
            for replacement in map(sp.Poly.as_expr, possible_replacements):
                new_system = deepcopy(curr_system)
                new_symbol, new_symbol_dot = new_system.variables.create_symbol_with_derivative()
                new_system.replace_subexpression(replacement, new_symbol)
                new_system._replacement_equations.append(sp.Eq(new_symbol, replacement))
                new_system._equations.append(sp.Eq(new_symbol_dot, new_system._calculate_Lie_derivative(replacement)).expand())

                system_queue.put(new_system)
            counter += 1

    def _debug_system_print(self, level: Optional[str]):
        if level is None or level == 'silent':
            pass
        elif level == 'info':
            print('-' * 100)
            print(f"Equations added: {len(self._replacement_equations)}")
            print(f"Last replacement: {self._replacement_equations[-1]}")
        elif level == 'debug':
            print('-' * 100)
            print(f"Equations added: {len(self._replacement_equations)}")
            print(f"Last replacement: {self._replacement_equations[-1]}")
            print('Equations:')
            print(self.equations)

        else:
            raise ValueError("debug value must be 'silent', 'info' or debug'")

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))
