import sympy as sp
from sympy.core.function import AppliedUndef
import hashlib
from copy import deepcopy
from typing import Callable, Set, Optional, List
from .AST_walk import find_non_polynomial
from .util import *
from .printer import str_common


class Variable(sp.Symbol):
    pass


class Parameter(sp.Symbol):
    pass


class VariablesHolder:
    """
    Class that manages variable storage.
    """

    def __init__(self, variables: Iterable[sp.Symbol],
                 parameter_variables: Optional[Set[sp.Symbol]] = None,
                 input_variables: Optional[Set[sp.Symbol]] = None,
                 new_var_base_name="w_",
                 start_new_vars_with=0):
        if parameter_variables is None:
            parameter_variables = set()
        if input_variables is None:
            input_variables = set()

        self._free_variables = list(variables)
        self._parameter_variables = parameter_variables
        self._input_variables = input_variables
        self.state_variables = list(variables)
        self._generated_variables = list()
        self._base_name = new_var_base_name
        self._start_id = start_new_vars_with

    @property
    def free(self):
        return self._free_variables

    @property
    def parameter(self):
        return self._parameter_variables

    @property
    def input(self):
        return self._input_variables

    @property
    def base_var_name(self):
        return self._base_name

    @base_var_name.setter
    def base_var_name(self, value: str):
        self._base_name = value

    @property
    def start_new_vars_with(self):
        return self._start_id

    @start_new_vars_with.setter
    def start_new_vars_with(self, value):
        self._start_id = value

    def create_variable(self) -> sp.Symbol:
        """
        Creates a new variable and stores it within itself.

        Example:
            .. math:: y_1 = holder.create\_variable()

        :return: Created variable
        """
        new_index = len(self._generated_variables) + self._start_id
        new_variable = sp.Symbol(self._base_name + "{%d}" % new_index)

        self._free_variables.append(new_variable)
        self.state_variables.append(new_variable)
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
        # _variables = _symbols.difference(_parameter_vars).difference(_input_vars)
        # _variables = set(filter(lambda v: r'\dot' not in str(v), _variables))
        _variables = remove_repeating([symbol_from_derivative(eq.lhs) for eq in equations])
        self.variables = VariablesHolder(_variables, _parameter_vars, _input_vars)

        self._equations_poly_degrees = dict()

    @property
    def equations(self) -> List[sp.Eq]:
        return self._equations

    @property
    def substitution_equations(self):
        return self._substitution_equations

    def to_poly_equations(self, inputs_ord: dict, optimize_parameters: set):
        """System should be already polynomial. Otherwise, throws Exception. You can check it by `is_polynomial` method
        :param optimize_parameters:
        """
        assert self.is_polynomial()
        inputs_ord_sym = {sp.Symbol(str_common(k)): v for k, v in inputs_ord.items()}
        opt_par_sym = {sp.Symbol(str_common(v)) for v in optimize_parameters}
        d_inputs = generate_derivatives(inputs_ord_sym)
        # TODO: Make explicit names for the highest order derivatives instead of 0
        coef_field = sp.FractionField(sp.QQ, list(map(str, self.variables.parameter)))
        input_ring_vars = list(dict.fromkeys(sp.flatten(d_inputs) + list(opt_par_sym)))
        R = coef_field[self.variables.state_variables + input_ring_vars]
        equations = [R.from_sympy(eq.rhs) for eq in self.equations]
        for i, v in enumerate(inputs_ord_sym.keys()):
            for dv in [g for g in R.gens if str(v) + '\'' in str(g)]:
                equations.append(dv)
            equations.append(R.zero)
        inputs_to_exclude = [tuple(R.from_sympy(v[-1])) for v in d_inputs if v[0] not in opt_par_sym]
        return equations, inputs_to_exclude

    def subs_expression(self, old: sp.Expr, new: sp.Expr):
        """Replace 'old' expression with 'new' expression for each equation."""
        for i in range(len(self._equations)):
            self._equations[i] = self._equations[i].subs(old, new)

    def replace_equation(self, old: sp.Expr, new: sp.Expr):
        for i in range(len(self._equations)):
            self._equations[i] = self._equations[i].replace(old, new)

    def expand_equations(self):
        """Apply SymPy 'expand' function to each of equation."""
        for i in range(len(self._equations)):
            self._equations[i] = sp.expand(self._equations[i])

    def is_polynomial(self) -> bool:
        for eq in self._equations:
            if not eq.args[1].is_polynomial(*self.variables.free, *self.variables.input):
                return False
        return True

    def auxiliary_equation_type_choose(self, auxiliary_eq_type: str) -> Callable:
        if auxiliary_eq_type == 'differential':
            return self.differential_auxiliary_equation_add
        elif auxiliary_eq_type == 'algebraic':
            return self.algebraic_auxiliary_equation_add
        else:
            raise ValueError("auxiliary_eq_type must be 'algebraic' or 'differential'")

    def differential_auxiliary_equation_add(self, new_var: sp.Symbol, substitution: sp.Expr) -> None:
        """
        Add differential auxiliary equation, generated by substitution.

        Substitution equation:
            .. math:: y = f(x)
        Generated equation:
            .. math:: \dot y = \dot f(x)

        :param new_var: left part, y
        :param substitution: right part, f(x)
        """
        if is_noninteger_positive_exp(substitution):
            self._add_noninteger_positive_subs_equation(new_var, substitution)
        else:
            self.replace_equation(substitution, new_var)
            self._substitution_equations.append(sp.Eq(new_var, substitution))
            self._equations.append(sp.Eq(
                make_derivative_symbol(new_var),
                self._calculate_Lie_derivative(substitution)).expand())

    def _add_noninteger_positive_subs_equation(self, new_var, substitution):
        int_part = int(substitution.exp)
        rem_part = substitution.exp - int_part
        coef = rem_part - 1
        var = substitution.args[0]

        print(f"{new_var} = {var ** (rem_part - 1)}")
        var_diff_eq = list(filter(lambda eq: eq.args[0] == make_derivative_symbol(var), self._equations))[0].rhs
        self._equations.append(sp.Eq(
            make_derivative_symbol(new_var), coef * new_var ** 2 * var_diff_eq / substitution * var ** int_part)
        )
        # x^2.3 -> z x^3
        self.replace_equation(substitution, new_var * var ** (int_part + 1))

    def algebraic_auxiliary_equation_add(self, new_variable: sp.Symbol, substitution: sp.Expr) -> None:
        """
        Add algebraic auxiliary equation, generated by substitution.

        substitution and added equation:
            .. math:: y = f(x)

        :param new_variable: left part, y
        :param substitution: right part, f(x)
        """
        self.subs_expression(substitution, new_variable)
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
            self.variables.input.add(input_var_dot)
            result += expr.diff(input_var) * input_var_dot
        return self._apply_substitutions(self._apply_substitutions(result).expand())

    def _apply_substitutions(self, expr: sp.Expr) -> sp.Expr:
        for left, right in map(lambda eq: eq.args, self._substitution_equations):
            expr = expr.subs(right, left)
        return expr

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
            print(rf"{eq.args[0]} = {sp.collect(eq.args[1], self.variables.free)}")

    def print_substitution_equations(self):
        for eq in self._substitution_equations:
            print(f"{eq.lhs} = {eq.rhs}")

    def __len__(self):
        return len(self._equations)

    def __repr__(self):
        return '\n'.join(map(lambda e: e.__repr__(), self._equations))

    def __str__(self):
        return '\n'.join(map(lambda e: e.__str__(), self._equations))


def polynomialize(system: Union[EquationSystem, List[Tuple[sp.Symbol, sp.Expr]]],
                  mode='differential', new_var_name="w_", start_new_vars_with=0) -> EquationSystem:
    """
    Transforms the system into polynomial form using variable substitution techniques.

    :param system: non-linear ODEs system
    :param mode: auxiliary equation form.
    :param new_var_name: base name for new variables. Example: new_var_name='w' => w0, w1, w2, ...
    :param start_new_vars_with: Initial index for new variables. Example: start_new_vars_with=3 => w3, w4, ...

    Mode
    -----------------
    **algebraic**
        adds auxiliary equations in form y = f(x, y)
    **differential**
         adds auxiliary equations in form y' = f(x, y)

    """
    if not isinstance(system, EquationSystem):
        system = eq_list_to_eq_system(system)
    system.variables.base_var_name = new_var_name
    system.variables.start_new_vars_with = start_new_vars_with
    if mode == 'algebraic':
        return _polynomialize_algebraic(system)
    elif mode == 'differential':
        return _polynomialize_differential(system)
    else:
        raise ValueError("mode must be 'algebraic' or 'differential")


def eq_list_to_eq_system(system: List[Tuple[sp.Symbol, sp.Expr]]) -> EquationSystem:
    lhs, rhs = zip(*system)
    params = set(reduce(lambda l, r: l | r, [eq.atoms(Parameter) for eq in rhs]))
    funcs = set(reduce(lambda l, r: l | r, [eq.atoms(AppliedUndef) for eq in rhs]))
    lhs_args = set(sp.flatten([eq.args for eq in lhs if not isinstance(eq, sp.Derivative)]))
    inputs = set(filter(lambda f: (f not in lhs) and (f not in lhs_args), funcs))
    spatial = funcs.difference(inputs)

    degrade_to_symbol = {s: sp.Symbol(str_common(s)) for s in funcs | inputs | params | set(lhs)}

    params_sym = [p.subs(degrade_to_symbol) for p in params]
    funcs_sym = [f.subs(degrade_to_symbol) for f in funcs]
    inputs_sym = [i.subs(degrade_to_symbol) for i in inputs]
    spacial_sym = [s.subs(degrade_to_symbol) for s in spatial]

    ders = derivatives([s.subs(degrade_to_symbol) for s in lhs])
    rhs_sym = [eq.subs(degrade_to_symbol) for eq in rhs]

    der_inputs = set(reduce(lambda l, r: l | r, [eq.atoms(sp.Derivative) for eq in rhs]))
    degrade_to_symbol.update({i: sp.Symbol(str_common(i)) for i in der_inputs})
    rhs_sym = [eq.subs(degrade_to_symbol) for eq in rhs]
    system = EquationSystem([sp.Eq(dx, fx) for dx, fx in zip(ders, rhs_sym)], params_sym, inputs_sym)
    system.variables._free_variables = list(set(system.variables.free + spacial_sym))
    return system


def _polynomialize_algebraic(system: EquationSystem) -> EquationSystem:
    result_system = deepcopy(system)
    while not result_system.is_polynomial():
        _polynomialize_algebraic_iter(result_system)

    return result_system


def _polynomialize_algebraic_iter(system: EquationSystem):
    for eq in system.equations:
        non_poly_elem = find_non_polynomial(eq.args[1], {})
        if non_poly_elem:
            new_variable = system.variables.create_variable()
            system.algebraic_auxiliary_equation_add(new_variable, non_poly_elem)
            break


def _polynomialize_differential(system: EquationSystem) -> EquationSystem:
    result_system = deepcopy(system)
    while not result_system.is_polynomial():
        _polynomialize_differential_iter(result_system)

    return result_system


def _polynomialize_differential_iter(system: EquationSystem):
    for eq in system.equations:
        non_poly_elem = find_non_polynomial(eq.args[1], set(system.variables.free) | system.variables.input)
        if non_poly_elem:
            new_variable = system.variables.create_variable()
            system.differential_auxiliary_equation_add(new_variable, non_poly_elem)
            break


def is_noninteger_positive_exp(expr: sp.Expr):
    return expr.is_Pow and expr.exp > 0 and (not expr.exp.is_integer)
