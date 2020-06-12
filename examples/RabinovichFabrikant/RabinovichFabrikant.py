import sympy as sp
import numpy as np
from qbee import *

sp.init_printing()

x, y, z = sp.symbols('x, y, z')
x_dot, y_dot, z_dot = derivatives('x, y, z')

a, b = sp.symbols('a, b')

system = EquationSystem([
    sp.Eq(x_dot, y * (z - 1 + x ** 2) + x),
    sp.Eq(y_dot, x * (3 * z + 1 - x ** 2) + y),
    sp.Eq(z_dot, -2 * z * (1 + x * y))
])


def ql_bfs_stat(_) -> int:
    try:
        ql_result = quadratic_linearize(system, method_optimal="bfs", limit_depth=3)
        return ql_result.statistics.steps
    except RuntimeError:
        return np.nan


def ql_iddfs_random(_) -> int:
    try:
        ql_result = quadratic_linearize(system, method_optimal="iddfs", heuristics='random', initial_max_depth=3, limit_depth=3)
        return ql_result.statistics.steps
    except RuntimeError:
        return np.nan


def ql_iddfs_frequent_first(_) -> int:
    try:
        ql_result = quadratic_linearize(system, method_optimal="iddfs", heuristics='frequent-first', initial_max_depth=3, limit_depth=3)
        return ql_result.statistics.steps
    except RuntimeError:
        return np.nan


def ql_iddfs_free_variables_count(_) -> int:
    try:
        ql_result = quadratic_linearize(system, method_optimal="iddfs", heuristics='free-variables-count', initial_max_depth=3, limit_depth=3)
        return ql_result.statistics.steps
    except RuntimeError:
        return np.nan


def ql_iddfs_auxiliary_equation_degree(_) -> int:
    try:
        ql_result = quadratic_linearize(system, method_optimal="iddfs", heuristics='auxiliary-equation-degree', initial_max_depth=3, limit_depth=3)
        return ql_result.statistics.steps
    except RuntimeError:
        return np.nan


def ql_iddfs_auxiliary_equation_ql_discrepancy(_) -> int:
    try:
        ql_result = quadratic_linearize(system, method_optimal="iddfs", heuristics='auxiliary-equation-ql-discrepancy', initial_max_depth=3,
                                        limit_depth=3)
        return ql_result.statistics.steps
    except RuntimeError:
        return np.nan


def ql_iddfs_summary_monomial_degree(_) -> int:
    try:
        ql_result = quadratic_linearize(system, method_optimal="iddfs", heuristics='summary-monomial-degree', initial_max_depth=3, limit_depth=3)
        return ql_result.statistics.steps
    except RuntimeError:
        return np.nan
