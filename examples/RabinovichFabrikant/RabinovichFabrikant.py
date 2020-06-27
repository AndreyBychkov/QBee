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


def quad_bfs_stat(_) -> int:
    try:
        quad_result = quadratize(system, method_optimal="bfs", limit_depth=3)
        return quad_result.statistics.steps
    except RuntimeError:
        return np.nan


def quad_iddls_random(_) -> int:
    try:
        quad_result = quadratize(system, method_optimal="iddls", heuristics='random', initial_max_depth=3, limit_depth=3)
        return quad_result.statistics.steps
    except RuntimeError:
        return np.nan


def quad_iddls_frequent_first(_) -> int:
    try:
        quad_result = quadratize(system, method_optimal="iddls", heuristics='frequent-first', initial_max_depth=3, limit_depth=3)
        return quad_result.statistics.steps
    except RuntimeError:
        return np.nan


def quad_iddls_free_variables_count(_) -> int:
    try:
        quad_result = quadratize(system, method_optimal="iddls", heuristics='free-variables-count', initial_max_depth=3, limit_depth=3)
        return quad_result.statistics.steps
    except RuntimeError:
        return np.nan


def quad_iddls_auxiliary_equation_degree(_) -> int:
    try:
        quad_result = quadratize(system, method_optimal="iddls", heuristics='auxiliary-equation-degree', initial_max_depth=3, limit_depth=3)
        return quad_result.statistics.steps
    except RuntimeError:
        return np.nan


def quad_iddls_auxiliary_equation_quadratic_discrepancy(_) -> int:
    try:
        quad_result = quadratize(system, method_optimal="iddls", heuristics='auxiliary-equation-quadratic-discrepancy', initial_max_depth=3,
                               limit_depth=3)
        return quad_result.statistics.steps
    except RuntimeError:
        return np.nan


def quad_iddls_summary_monomial_degree(_) -> int:
    try:
        quad_result = quadratize(system, method_optimal="iddls", heuristics='summary-monomial-degree', initial_max_depth=3, limit_depth=3)
        return quad_result.statistics.steps
    except RuntimeError:
        return np.nan
