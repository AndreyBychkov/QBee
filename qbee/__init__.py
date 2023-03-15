"""
QBee - Package for transforming Systems of Differential Equations to quadratic form.
"""

__version__ = "0.8.0"

from .quadratization import polynomialize_and_quadratize, quadratize
from .polynomialization import polynomialize
from .dimension_agnostic import quadratize_dimension_agnostic
from .quadratization import pruning_by_nodes_processed, pruning_by_vars_number, pruning_by_declining_variables, \
    pruning_by_elapsed_time, pruning_by_nodes_without_quadratization_found, \
    default_pruning_rules, without_variables
from .selection import default_generation, generation_semidiscretized, default_scoring, smd_scoring, aeqd_scoring
from .printer import print_qbee, str_qbee
from .util import functions, parameters, INDEPENDENT_VARIABLE
