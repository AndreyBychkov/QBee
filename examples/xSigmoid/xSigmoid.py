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
    quad_result = quadratize(system, method_optimal="bfs", limit_depth=2)
    return quad_result.statistics.steps


def quad_iddfs_random(_) -> int:
    quad_result = quadratize(system, method_optimal="iddfs", heuristics='random', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddfs_frequent_first(_) -> int:
    quad_result = quadratize(system, method_optimal="iddfs", heuristics='frequent-first', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddfs_free_variables_count(_) -> int:
    quad_result = quadratize(system, method_optimal="iddfs", heuristics='free-variables-count', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddfs_auxiliary_equation_degree(_) -> int:
    quad_result = quadratize(system, method_optimal="iddfs", heuristics='auxiliary-equation-degree', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps


def quad_iddfs_auxiliary_equation_quadratic_discrepancy(_) -> int:
    quad_result = quadratize(system, method_optimal="iddfs", heuristics='auxiliary-equation-quadratic-discrepancy', initial_max_depth=2,
                           limit_depth=2)
    return quad_result.statistics.steps


def quad_iddfs_summary_monomial_degree(_) -> int:
    quad_result = quadratize(system, method_optimal="iddfs", heuristics='summary-monomial-degree', initial_max_depth=2, limit_depth=2)
    return quad_result.statistics.steps
