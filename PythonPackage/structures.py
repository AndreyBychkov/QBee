import sympy as sp
import hashlib

from functools import reduce
from typing import List, Iterable, Optional, Callable, Tuple, Set
from util import polynomial_replace, get_possible_replacements
from statistics import *


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


class SymbolsHolder:
    """
    Class that manages symbol storage.
    """

    def __init__(self, symbols: Iterable):
        self._original_symbols = list(symbols)
        self._added_symbols = list()
        self._base_name = "y_"

    def create_symbol(self) -> sp.Symbol:
        """
        Creates a new symbol and stores it within itself.

        Example:
            .. math:: y_1 = holder.create\_symbol()

        :return: Created symbol
        """
        if not self._added_symbols:
            new_symbol = sp.Symbol(self._base_name + "{0}")
        else:
            last = str(self.get_last_added())
            new_index = int(last.split("_")[1][1:-1]) + 1
            new_symbol = sp.Symbol(self._base_name + "{%d}" % new_index)

        self._added_symbols.append(new_symbol)
        return new_symbol

    def create_symbol_with_derivative(self) -> Tuple[sp.Symbol, sp.Symbol]:
        """
        Creates new a symbol with its derivative and stores them within itself.

        Example:
            .. math:: y_1, \dot{y}_1 = holder.create\_symbol\_with\_derivative()

        :returns: Created symbol with derivative
        """
        new_symbol = self.create_symbol()
        new_symbol_der = make_derivative_symbol(new_symbol)
        return new_symbol, new_symbol_der

    def add_symbols(self, symbols: Iterable):
        for new_symbol in symbols:
            self._added_symbols.append(new_symbol)

    def get_last_added(self) -> sp.Symbol:
        return self._added_symbols[-1]

    def get_symbols(self) -> List[sp.Symbol]:
        return self._original_symbols + self._added_symbols


class EquationSystem:
    def __init__(self, equations: List[sp.Eq],
                 parameter_variables: Iterable[sp.Symbol] = None,
                 input_variables: Iterable[sp.Symbol] = None):
        self._equations = equations.copy()
        self._original_equation_indexes = list(range(len(equations)))
        self._replacement_equations = list()

        self._equations_poly_degrees = dict()
        if self.is_polynomial():
            self._compute_equations_poly_degrees()

        self.variables = SymbolsHolder(reduce(set.union, map(lambda e: e.free_symbols, equations)))
        self._parameter_vars = set(parameter_variables) if parameter_variables is not None else set()
        self._input_vars = set(input_variables) if input_variables is not None else set()

        self.statistics: EvaluationStatistics = EvaluationStatistics(0, 0, "empty")

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

    def replace_expression(self, old: sp.Expr, new: sp.Expr):
        """Replace 'old' expression with 'new' expression for each equation."""
        for i in range(len(self._equations)):
            self._equations[i] = self._equations[i].subs(old, new)

    def replace_monomial(self, old: sp.Expr, new: sp.Expr):
        """If any expression in system is divisible on 'old', replace it with 'new'"""
        for i, eq in enumerate(self._equations):
            self._equations[i] = sp.Eq(eq.args[0], polynomial_replace(eq.args[1], old, new))

    def expand_equations(self):
        """Apply SymPy 'expand' function to each of equation."""
        for i in range(len(self._equations)):
            self._equations[i] = sp.expand(self._equations[i])

    def update_poly_degrees(self):
        """Update system variable _equations_poly_degrees"""
        self._compute_equations_poly_degrees()

    def get_possible_replacements(self, count_sorted=False) -> Tuple[sp.Poly]:
        right_equations = self._get_right_equations()
        return get_possible_replacements(right_equations, self._parameter_vars, count_sorted=count_sorted)

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
            if sp.total_degree(mon, *poly.free_symbols.difference(self._parameter_vars)) > 2:
                return False
        return True

    def auxiliary_equation_type_choose(self, auxiliary_eq_type: str) -> Callable:
        if auxiliary_eq_type == 'differential':
            return self.differential_auxiliary_equation_add
        elif auxiliary_eq_type == 'algebraic':
            return self.algebraic_auxiliary_equation_add
        else:
            raise ValueError("auxiliary_eq_type must be 'algebraic' or 'differential'")

    def differential_auxiliary_equation_add(self, new_symbol: sp.Symbol, replacement: sp.Expr, is_polynomial_replacement=True) -> None:
        """
        Add differential auxiliary equation, generated by replacement.

        Replacement equation:
            .. math:: y = f(x)
        Generated equation:
            .. math:: \dot y = \dot f(x)

        :param new_symbol: left part, y
        :param replacement: right part, f(x)
        :param is_polynomial_replacement: is replacement is polynomial
        """
        new_symbol_dot = make_derivative_symbol(new_symbol)
        if is_polynomial_replacement:
            self.replace_monomial(replacement, new_symbol)
        else:
            self.replace_expression(replacement, new_symbol)
        self._replacement_equations.append(sp.Eq(new_symbol, replacement))
        self._equations.append(sp.Eq(new_symbol_dot, self._calculate_Lie_derivative(replacement)).expand())

    def algebraic_auxiliary_equation_add(self, new_symbol: sp.Symbol, replacement: sp.Expr, is_polynomial_replacement=True) -> None:
        """
        Add algebraic auxiliary equation, generated by replacement.

        Replacement and added equation:
            .. math:: y = f(x)

        :param new_symbol: left part, y
        :param replacement: right part, f(x)
        :param is_polynomial_replacement: is replacement is polynomial
        """
        if is_polynomial_replacement:
            self.replace_monomial(replacement, new_symbol)
        else:
            self.replace_expression(replacement, new_symbol)
        self._replacement_equations.append(sp.Eq(new_symbol, replacement))
        self._equations.append(sp.Eq(new_symbol, replacement).expand())

    def _calculate_Lie_derivative(self, expr: sp.Expr) -> sp.Expr:
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

    def _replace_from_replacement_equations(self, expr: sp.Expr) -> sp.Expr:
        new_expr = expr.copy()
        for left, right in map(lambda eq: eq.args, self._replacement_equations):
            new_expr = new_expr.subs(right, left)
        return new_expr

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

    def _get_right_equations(self):
        return list(map(lambda eq: eq.args[1], self._equations))

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))


class QuadraticLinearizationResult:
    def __init__(self, system: EquationSystem, statistics: EvaluationStatistics, replacements: Set[sp.Eq]):
        self.system = system
        self.statistics = statistics
        self.replacements = replacements
