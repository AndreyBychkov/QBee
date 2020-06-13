import sympy as sp
from qbee import *

sp.init_printing()

x = sp.symbols('x')
dot_x = derivatives(x)

system = EquationSystem([
    sp.Eq(dot_x, x / (1 + sp.exp(x)))
])

system = polynomialize(system)


def ql_bfs_stat(_) -> int:
    ql_result = quadratize(system, method_optimal="bfs", limit_depth=2)
    return ql_result.statistics.steps


def ql_iddfs_random(_) -> int:
    ql_result = quadratize(system, method_optimal="iddfs", heuristics='random', initial_max_depth=2, limit_depth=2)
    return ql_result.statistics.steps


def ql_iddfs_frequent_first(_) -> int:
    ql_result = quadratize(system, method_optimal="iddfs", heuristics='frequent-first', initial_max_depth=2, limit_depth=2)
    return ql_result.statistics.steps


def ql_iddfs_free_variables_count(_) -> int:
    ql_result = quadratize(system, method_optimal="iddfs", heuristics='free-variables-count', initial_max_depth=2, limit_depth=2)
    return ql_result.statistics.steps


def ql_iddfs_auxiliary_equation_degree(_) -> int:
    ql_result = quadratize(system, method_optimal="iddfs", heuristics='auxiliary-equation-degree', initial_max_depth=2, limit_depth=2)
    return ql_result.statistics.steps


def ql_iddfs_auxiliary_equation_ql_discrepancy(_) -> int:
    ql_result = quadratize(system, method_optimal="iddfs", heuristics='auxiliary-equation-ql-discrepancy', initial_max_depth=2,
                           limit_depth=2)
    return ql_result.statistics.steps


def ql_iddfs_summary_monomial_degree(_) -> int:
    ql_result = quadratize(system, method_optimal="iddfs", heuristics='summary-monomial-degree', initial_max_depth=2, limit_depth=2)
    return ql_result.statistics.steps
