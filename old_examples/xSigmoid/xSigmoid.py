import sympy as sp
from qbee import *

sp.init_printing()

x = sp.symbols('x')
dot_x = derivatives(x)

system = EquationSystem([
    sp.Eq(dot_x, x / (1 + sp.exp(x)))
])

system = polynomialize(system)


def quad_bfs_stat(_) -> int:
    quad_result = quadratize(system, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddls_random(_) -> int:
    quad_result = quadratize(system, heuristics='random', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddls_frequent_first(_) -> int:
    quad_result = quadratize(system, heuristics='frequent-first', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddls_free_variables_count(_) -> int:
    quad_result = quadratize(system, heuristics='free-variables-count', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddls_auxiliary_equation_degree(_) -> int:
    quad_result = quadratize(system, heuristics='auxiliary-equation-degree', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddls_auxiliary_equation_quadratic_discrepancy(_) -> int:
    quad_result = quadratize(system, heuristics='auxiliary-equation-quadratic-discrepancy', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddls_summary_monomial_degree(_) -> int:
    quad_result = quadratize(system, heuristics='summary-monomial-degree', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps
