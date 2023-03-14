"""
QBee - Package for transforming Systems of Differential Equations to quadratic form.
"""

__version__ = "0.8.0"

from .polynomialization import polynomialize
from .quadratization import quadratize, polynomialize_and_quadratize, PolynomialSystem, \
    pruning_by_nodes_processed, pruning_by_vars_number, pruning_by_quadratic_upper_bound, pruning_by_squarefree_graphs, \
    pruning_by_best_nvars, pruning_by_elapsed_time, pruning_by_nodes_without_quadratization_found, \
    default_pruning_rules, without_variables
from .dimension_agnostic import quadratize_dimension_agnostic
from .selection import *
from .printer import print_qbee, str_qbee
from .util import functions, parameters, multivariable_functions, Parameter, INDEPENDENT_VARIABLE
