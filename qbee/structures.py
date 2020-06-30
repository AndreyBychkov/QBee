import sympy as sp
import hashlib

from functools import reduce
from typing import List, Iterable, Optional, Callable, Tuple, Set
from .util import *


def derivatives(names) -> Tuple[sp.Symbol]:
    """
    Add dot to input symbols.

    :param names: input symbols. Can be represented in different ways.
    :return: Input symbols wrapper with dots.

    Example:
        .. math:: \dot{x}, \dot{y} = derivatives([x, y])

    Names variants
    -----------------
    **symbol**
        derivatives('x'), derivatives(x: Symbol)
    **string with delimiters**
        derivatives('x, y, z'), derivatives('x y z')
    **iterable of symbols**
        derivatives([x, y, z]), derivatives(['x', 'y', 'z'])
    """
    if not isinstance(names, Iterable) and isinstance(names, sp.Symbol):
        return (make_derivative_symbol(names))

    if isinstance(names, Iterable) and reduce(lambda a, b: a and b, map(lambda elem: isinstance(elem, sp.Symbol), names)):
        return tuple(map(make_derivative_symbol, names))

    symbols_output = sp.symbols(names)
    if isinstance(symbols_output, sp.Symbol):
        return make_derivative_symbol(symbols_output)
    else:
        return tuple(map(make_derivative_symbol, sp.symbols(names)))


def make_derivative_symbol(symbol) -> sp.Symbol:
    """
    Makes symbol with the dot from input symbol.

    If symbols with index must be in form 's_{index}'.

    Example:
        .. math:: \dot{x} = make\_derivative\_symbol(x)

    """
    str_symbol = str(symbol)
    if '_' in str_symbol:
        name, index, *other = str(symbol).split('_')
        return sp.Symbol(rf'\dot {name}' + '_' + '%s' % index)
    else:
        return sp.Symbol(rf'\dot {str_symbol}')


class VariablesHolder:
    """
    Class that manages variable storage.
    """

    def __init__(self, variables: Iterable[sp.Symbol],
                 parameter_variables: Optional[Set[sp.Symbol]] = None,
                 input_variables: Optional[Set[sp.Symbol]] = None):
        if parameter_variables is None:
            parameter_variables = set()
        if input_variables is None:
            input_variables = set()

        self._free_variables = list(variables)
        self._parameter_variables = parameter_variables
        self._input_variables = input_variables
        self._generated_variables = list()
        self._base_name = "y_"

    @property
    def free(self):
        return self._free_variables

    @property
    def parameter(self):
        return self._parameter_variables

    @property
    def input(self):
        return self._input_variables

    def create_variable(self) -> sp.Symbol:
        """
        Creates a new variable and stores it within itself.

        Example:
            .. math:: y_1 = holder.create\_variable()

        :return: Created variable
        """
        new_index = len(self._generated_variables)
        new_variable = sp.Symbol(self._base_name + "{%d}" % new_index)

        self._free_variables.append(new_variable)
        self._generated_variables.append(new_variable)
        return new_variable

    def create_variable_with_derivative(self) -> Tuple[sp.Symbol, sp.Symbol]:
        """
        Creates new a variable with its derivative and stores them within itself.

        Example:
            .. math:: y_1, \dot{y}_1 = holder.create\_variable\_with\_derivative()

        :returns: Created variable with derivative
        """
        new_variable = self.create_variable()
        new_variable_der = make_derivative_symbol(new_variable)
        return new_variable, new_variable_der


class EquationSystem:
    def __init__(self, equations: List[sp.Eq],
                 parameter_variables: Iterable[sp.Symbol] = None,
                 input_variables: Iterable[sp.Symbol] = None):
        self._equations = equations.copy()
        self._original_equation_indexes = list(range(len(equations)))
        self._substitution_equations = list()

        self.expand_equations()

        _symbols = reduce(set.union, map(lambda e: e.free_symbols, equations))
        _parameter_vars = set(parameter_variables) if parameter_variables is not None else set()
        _input_vars = set(input_variables) if input_variables is not None else set()
        _variables = _symbols.difference(_parameter_vars).difference(_input_vars)
        _variables = set(filter(lambda v: r'\dot' not in str(v), _variables))
        self.variables = VariablesHolder(_variables, _parameter_vars, _input_vars)

        self._equations_poly_degrees = dict()
        if self.is_polynomial():
            self._compute_equations_poly_degrees()

    @property
    def equations(self) -> List[sp.Eq]:
        return self._equations

    @property
    def equations_hash(self) -> bytes:
        return hashlib.md5(str(self._equations).encode('utf-8')).digest()

    @property
    def monomials(self) -> Tuple[sp.Expr]:
        """Sequential non-unique monomials of system"""
        assert self.is_polynomial("full")

        return reduce(lambda a, b: a + b, map(sp.Add.make_args, self._get_right_equations()))

    @property
    def substitution_equations(self):
        return self._substitution_equations

    def replace_expression(self, old: sp.Expr, new: sp.Expr):
        """Replace 'old' expression with 'new' expression for each equation."""
        for i in range(len(self._equations)):
            self._equations[i] = self._equations[i].subs(old, new)

    def replace_monomial(self, old: sp.Expr, new: sp.Expr):
        """If any expression in system is divisible on 'old', replace it with 'new'"""
        for i, eq in enumerate(self._equations):
            self._equations[i] = sp.Eq(eq.args[0], polynomial_subs(eq.args[1], old, new))

    def expand_equations(self):
        """Apply SymPy 'expand' function to each of equation."""
        for i in range(len(self._equations)):
            self._equations[i] = sp.expand(self._equations[i])

    def update_poly_degrees(self):
        """Update system variable _equations_poly_degrees"""
        self._compute_equations_poly_degrees()

    def get_possible_substitutions(self, count_sorted=False) -> Tuple[sp.Poly]:
        right_equations = self._get_right_equations()
        return get_possible_substitutions(right_equations, set(self.variables.free), count_sorted=count_sorted)

    def get_monomial_decompositions(self):
        right_equations = self._get_right_equations()
        return get_monomial_decompositions(right_equations, set(self.variables.free))

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
            if not self.equations[i].args[1].is_polynomial(*self.variables.free):
                return False
        return True

    def _is_polynomial_full(self) -> bool:
        for eq in self._equations:
            if not eq.args[1].is_polynomial(*self.variables.free):
                return False
        return True

    def is_quadratic(self, mode='full') -> bool:
        if mode == 'original':
            return self._is_quadratic_original()
        elif mode == 'full':
            return self._is_quadratic_full()
        else:
            raise ValueError("mode must be 'original' or 'full'.")

    def _is_quadratic_full(self) -> bool:
        for eq in self._equations:
            if not self._is_poly_quadratic(eq.args[1]):
                return False
        return True

    def _is_quadratic_original(self) -> bool:
        for i in self._original_equation_indexes:
            if not self._is_poly_quadratic(self._equations[i].args[1]):
                return False
        return True

    def _is_poly_quadratic(self, poly: sp.Expr) -> bool:
        monomials = sp.Add.make_args(poly)
        for mon in monomials:
            if sp.total_degree(mon, *poly.free_symbols.difference(self.variables.parameter)) > 2:
                return False
        return True

    def auxiliary_equation_type_choose(self, auxiliary_eq_type: str) -> Callable:
        if auxiliary_eq_type == 'differential':
            return self.differential_auxiliary_equation_add
        elif auxiliary_eq_type == 'algebraic':
            return self.algebraic_auxiliary_equation_add
        else:
            raise ValueError("auxiliary_eq_type must be 'algebraic' or 'differential'")

    def differential_auxiliary_equation_add(self, new_variable: sp.Symbol, substitution: sp.Expr, is_polynomial_substitution=True) -> None:
        """
        Add differential auxiliary equation, generated by substitution.

        Substitution equation:
            .. math:: y = f(x)
        Generated equation:
            .. math:: \dot y = \dot f(x)

        :param new_variable: left part, y
        :param substitution: right part, f(x)
        :param is_polynomial_substitution: is substitution is polynomial
        """
        new_variable_dot = make_derivative_symbol(new_variable)
        if is_polynomial_substitution:
            self.replace_monomial(substitution, new_variable)
        else:
            self.replace_expression(substitution, new_variable)
        self._substitution_equations.append(sp.Eq(new_variable, substitution))
        self._equations.append(sp.Eq(new_variable_dot, self._calculate_Lie_derivative(substitution)).expand())

    def algebraic_auxiliary_equation_add(self, new_variable: sp.Symbol, substitution: sp.Expr, is_polynomial_substitution=True) -> None:
        """
        Add algebraic auxiliary equation, generated by substitution.

        substitution and added equation:
            .. math:: y = f(x)

        :param new_variable: left part, y
        :param substitution: right part, f(x)
        :param is_polynomial_substitution: is substitution is polynomial
        """
        if is_polynomial_substitution:
            self.replace_monomial(substitution, new_variable)
        else:
            self.replace_expression(substitution, new_variable)
        self._substitution_equations.append(sp.Eq(new_variable, substitution))
        self._equations.append(sp.Eq(new_variable, substitution).expand())

    def _calculate_Lie_derivative(self, expr: sp.Expr) -> sp.Expr:
        """Calculates Lie derivative using chain rule."""
        result = sp.Integer(0)
        for var in expr.free_symbols.difference(self.variables.parameter).difference(self.variables.input):
            var_diff_eq = list(filter(lambda eq: eq.args[0] == make_derivative_symbol(var), self._equations))[0]
            var_diff = var_diff_eq.args[1]
            result += expr.diff(var) * var_diff
        for input_var in expr.free_symbols.intersection(self.variables.input):
            input_var_dot = make_derivative_symbol(input_var)
            result += expr.diff(input_var) * input_var_dot
        return self._apply_substitutions(self._apply_substitutions(result).expand())

    def _apply_substitutions(self, expr: sp.Expr) -> sp.Expr:
        for left, right in map(lambda eq: eq.args, self._substitution_equations):
            expr = expr.subs(right, left)
        return expr

    def _compute_equations_poly_degrees(self):
        for var_dot, right_expr in map(lambda eq: eq.args, self._equations):
            if right_expr.is_Number:
                self._equations_poly_degrees[var_dot] = 0
            else:
                self._equations_poly_degrees[var_dot] = sp.Poly(right_expr).total_degree()

    def _debug_system_print(self, level: Optional[str]) -> None:
        if level is None or level == 'silent':
            pass
        elif level == 'info':
            print('-' * 100)
            print(f"Equations added: {len(self._substitution_equations)}")
            print(f"Last substitution: {self._substitution_equations[-1]}")
        elif level == 'debug':
            print('-' * 100)
            print(f"Equations added: {len(self._substitution_equations)}")
            print(f"Last substitution: {self._substitution_equations[-1]}")
            print('Equations:')
            print(self.equations)

        else:
            raise ValueError("debug value must be 'silent', 'info' or debug'")

    def _get_right_equations(self):
        return list(map(lambda eq: eq.args[1], self._equations))

    def _get_right_substitution_equations(self):
        return list(map(lambda eq: eq.args[1], self._substitution_equations))

    def print(self, mode: str = 'simple'):
        if mode == "simple":
            self._print_simple()
        elif mode == 'latex':
            self._print_latex()
        elif mode == 'sympy':
            print(str(self))
        else:
            raise AttributeError(f"mode {mode} is not valid. Use correct mode.")

    def _print_latex(self):
        print(r'\begin{array}{ll}')
        for eq in self._equations:
            print('\t' + rf"{eq.args[0]} = {sp.latex(sp.collect(eq.args[1], self.variables.free))}" + r'\\')
        print(r'\end{array}')

    def _print_simple(self):
        for eq in self.equations:
            print(rf"{symbol_from_derivative(eq.args[0])}' = {sp.collect(eq.args[1], self.variables.free)}")

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))


class EvaluationStatistics:
    """Statistics structure for quadratization"""

    def __init__(self, steps: int, depth: int, method_name: str):
        """
        :param steps: number of computed systems
        :param depth: number of substitutions to make system quadratic
        :param method_name: name of the algorithm
        """
        self.steps = steps
        self.method_name = method_name
        self.depth = depth

    def __repr__(self):
        return '\n'.join([
            f"steps: {self.steps}",
            f"Method's name: {self.method_name}",
            f"depth: {self.depth}"
        ])


class QuadratizationResult:
    def __init__(self, system: EquationSystem, statistics: EvaluationStatistics, substitutions: Tuple[sp.Expr]):
        self.system = system
        self.statistics = statistics
        self.substitutions = substitutions
