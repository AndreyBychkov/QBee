"""
QBee - Package for transforming Systems of Differential Equations to quadratic form.
"""

__version__ = "0.6.0"

import sympy as sp
from typing import Iterable
from .polynomialization import polynomialize, EquationSystem, Parameter
from .quadratization import quadratize, polynomialize_and_quadratize, BranchAndBound, PolynomialSystem, Pruning, \
    pruning_by_nodes_processed, pruning_by_vars_number, pruning_by_quadratic_upper_bound, pruning_by_squarefree_graphs, \
    pruning_by_best_nvars, pruning_by_elapsed_time, pruning_by_nodes_without_quadratization_found, \
    default_pruning_rules, \
    without_variables
from .dimension_agnostic import quadratize_dimension_agnostic
from .selection import *
from .printer import print_qbee

INDEPENDENT_VARIABLE = sp.Symbol("_t")


def functions(names, *, laurent=True, **kwargs):
    """Dependent and input variables in differential equations. The syntax is identical to `sympy.symbols`_.

    Examples
        >>> x, y, z = functions('x, y, z')

    .. _sympy.symbols: https://docs.sympy.org/latest/modules/core.html?highlight=symbols#sympy.core.symbol.symbols
    """
    t = INDEPENDENT_VARIABLE
    funcs = sp.symbols(names, cls=sp.Function, is_laurent=laurent, **kwargs)
    if isinstance(funcs, Iterable):
        return [f(t) for f in funcs]
    return funcs(t)


def multivariable_functions(names, args, *, laurent=True, **kwargs):
    funcs = sp.symbols(names, cls=sp.Function, is_laurent=laurent, **kwargs)
    if isinstance(funcs, Iterable):
        return [f(*args) for f in funcs]
    return funcs(*args)


def parameters(names, **kwargs):
    """Parameter variables in differential equations. The syntax is identical to `sympy.symbols`_.

        Examples
            >>> p, k = parameters('p, k')

        .. _sympy.symbols: https://docs.sympy.org/latest/modules/core.html?highlight=symbols#sympy.core.symbol.symbols
        """
    return sp.symbols(names, cls=Parameter, **kwargs)
