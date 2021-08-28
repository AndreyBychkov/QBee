import sympy as sp
from .polynomialization import polynomialize, EquationSystem, Variable, Parameter
from .quadratization import BranchAndBound, PolynomialSystem, Pruning, \
    pruning_by_nodes_processed, pruning_by_vars_number, pruning_by_quadratic_upper_bound, pruning_by_squarefree_graphs, \
    pruning_by_best_nvars, pruning_by_elapsed_time, pruning_by_excluding_variables, default_pruning_rules, \
    quadratize, polynomialize_and_quadratize
from .selection import *
from .util import derivatives


def functions(names, **kwargs):
    return sp.symbols(names, cls=Variable, **kwargs)


def parameters(names, **kwargs):
    return sp.symbols(names, cls=Parameter, **kwargs)
