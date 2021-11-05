"""
QBee - Package for transforming Systems of Differential Equations to quadratic form.
"""

__version__ = "0.6.0"

import sympy as sp
from .polynomialization import polynomialize, EquationSystem, Variable, Parameter
from .quadratization import BranchAndBound, PolynomialSystem, Pruning, \
    pruning_by_nodes_processed, pruning_by_vars_number, pruning_by_quadratic_upper_bound, pruning_by_squarefree_graphs, \
    pruning_by_best_nvars, pruning_by_elapsed_time, without_variables, default_pruning_rules, \
    quadratize, polynomialize_and_quadratize
from .selection import *
from .util import derivatives


def functions(names, **kwargs):
    """Dependent and input variables in differential equations. The syntax is identical to `sympy.symbols`_.

    Examples
        >>> x, y, z = functions('x, y, z')

    .. _sympy.symbols: https://docs.sympy.org/latest/modules/core.html?highlight=symbols#sympy.core.symbol.symbols
    """
    return sp.symbols(names, cls=Variable, **kwargs)


def parameters(names, **kwargs):
    """Parameter variables in differential equations. The syntax is identical to `sympy.symbols`_.

        Examples
            >>> p, k = parameters('p, k')

        .. _sympy.symbols: https://docs.sympy.org/latest/modules/core.html?highlight=symbols#sympy.core.symbol.symbols
        """
    return sp.symbols(names, cls=Parameter, **kwargs)
